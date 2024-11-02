#!/usr/bin/env python3

import io
import os
import sys
import glob
import socket
import selectors
import numpy


# version 0.1       : start, use fewer modules & codes as possible as we can
# version 0.2       : receive command line inputs
# version 0.3       : BC, remove selectors.EVENT_WRITE
# version 0.4       : check receives
# version 0.5       : integration, avoid intermediate database file
# version 0.6       : integrate numpyfy_dict to qmtxt2db
# version 0.6.1     : get PORT from environment variable, if not, exit
# version 0.7       : multiple models
# version 0.8       : BC, integrate aenet
# version 0.8.1     : for aenet server, add choice for unit conversions
# version 0.8.2     : elevate the numpy package and cleanup the formats


"""ML server to speed up schnetpack & aenet prediction

this server assumes model file and database file already exist under
current path, server is initialized once but makes prediction many times
during its live time whenever it receives 'spk' or 'aenet' bytes.

To be safe, this server is built up only on internal network (works only
inside machine), and it can only listen on two ports at the same time.

Although it can send and receive messages, only 'spk' or 'aenet' make sense,
otherwise, it will terminate itself immediately, and print out this wrong bytes.

this server is an integration server, with many updates and revisions,
which can directly parse qm.txt file and make prediction on-the-fly during
its runtime, to avoid the generation of intermediate files.


logical:
    for schnetpack:
        qmtxt2db: in heavy duty, parse qm.txt file
        AtomsData: revised source codes
        AtomsLoader: revised source codes

    for aenet:
        qmtxt2atoms: convert qmfile to atoms list
        ANNCalculatorMP: multiprocessing python-aenet-ANN feeder


main methods(same name but different on functions):
    stdselection: selection based on standard deviation
    accept: build up incoming connections, and register them to selectors
    myconnect: receive messages and process accordingly
    predict: main, output results to new "ELEC.energy" file


note:
    this script has to be run directly, it cannot be imported by third module(s)
    coding in this way break PEP8 recommends
"""


USAGE = """myserver usage
    ./myserver  [spk|aenet]
"""

if len(sys.argv) != 2:
    print('Error: wrong using')
    print(USAGE)
    sys.exit(-1)
PROGRAM = sys.argv[1].lower()
if PROGRAM not in ['spk', 'aenet']:
    print('Error: wrong program: {:}'.format(sys.argv[1]))
    sys.exit(-1)
print(f'Note: program to choose: {PROGRAM}')


# to get it exact IP, execute `ping localhost'
# normally, it will be set to `127.0.0.1'
HOST = 'localhost'


# deprecated
#
# Ports defined:
#   0-1023      : the Well Known Ports, also referred to as System Ports.
#   1024-49151  : the Registered Ports, also known as User Ports.
#   49152-65535 : the Dynamic Ports, also referred to as the Private Ports.
#
# thus, please set this number in range 1024~65535
#
# now, this number is set at its toplevel script: `pmfcmd'
#PORT = 9999


# valid when program spk is chosen
#
# name format, any model starts with this string will be loaded
#
# please do not use white space on the name of model files
SPKMD = 'best_model'


# valid when program aenet is chosen
#
# the number of dicts in this list indicates the number of trained models
#
# format is in python-dict: `trained-atom-type : trained-model-file'
#
# for each its entry, additional key `unit' can be used,
# values choice from: `hartree' and `kcal', default is `kcal'.
# conversion is: 1 hartree = 627.509469 kcal/mol
#
# this is useful when trained model is using the different unit
AENETMDS = [
    {
        'unit': 'hartree',          # unit will be hartree
        'C': 'C.16tw-16tw-8tw.ann',
        'H': 'H.16tw-16tw-8tw.ann',
        'O': 'O.16tw-16tw-8tw.ann',
    },
]


# standard deviation, the selection range will be:
#   P(mean-SCALE*std <= X <= mean+SCALE*std)
SCALE = 1.2


# whether print more info, y/yes, (only 'y' or 'yes' works), else not.
# please take care of variable's format, string 'y' or 'yes' has to be quoted.
# due to server is running in background in the different process,
# thus if debug mode is on, a new file called "debuginfo.txt" will be generated,
# which will include the name of loaded models (defined by SPKMD or predict*in)
# and their corresponding predicted values
DEBUG = 'n'


# coordinates file that BOSS will generate
QMFILE = 'qm.txt'

# bytes to communicate
BYTES = 8

#
# END User edits
#

# patch for v0.61
PORT = os.getenv('PORT')
bo = False
if PORT is None:
    bo = True
elif not isinstance(PORT,int):
    try:
        PORT = int(PORT)
    except:
        bo = True
if bo:
    print('Warning: PORT is not set')
    print('       : script can only be called inside other script')
    print('       : now it is not a stand-alone server')
    sys.exit(-1)


sock = socket.socket()
rst = sock.connect_ex((HOST,PORT))
if rst == 0:
    print('Warning: port is already opened: {:}:{:}'.format(HOST,PORT))
    print('       : highly like previous run is not properly terminated')
    print('       : please execute < lsof -i4 | grep -i listen > to check')
    print('       : outputs will be in format: < python3 PID user 3u IPV4... >')
    print('       : get second entry (integer) process ID < PID >')
    print('Warning: Be extremely careful, make sure 100% confident')
    print('       : kill it by: < kill PID >')
    print('       : Or, change port number to other value bigger than 1023')
    print('       : at the same time please update CSHELL script and try again')
    sys.exit(-1)


# process DEBUG to True or False
DEBUG = True if str(DEBUG).lower() in ['y','yes'] else False


# since BOSS qm.txt file is using numbers indicate true atom type
# thus we separately defined three constants to do convension
PERODIC_TABLE = {
    1:'H', 5:'B', 6:'C', 7:'N', 8:'O', 9:'F', 15:'P', 16:'S', 17:'Cl',
    35:'Br', 53:'I',
}
PT2 = {}
ATOMTYPES = []
for k,v in PERODIC_TABLE.items():
    PT2[str(k)] = v
    ATOMTYPES.append(v)


def stdselection(enelist,scale=None):
    proenelist = [i for i in enelist if isinstance(i,(float,int))]
    if not proenelist: return "FAILED"
    if len(proenelist) <= 2: return sum(proenelist)/len(proenelist)
    e = numpy.mean(proenelist)
    if scale is None:
        return e
    s = numpy.std(proenelist)
    h = e + s*scale
    l = e - s*scale
    if h < l: h,l = l,h
    ls = [i for i in proenelist if i <= h and i >= l]
    if not ls:
        return e
    return sum(ls)/len(ls)


if PROGRAM == 'aenet':

    for md in AENETMDS:
        for k,v in md.items():
            if k != 'unit' and not os.path.isfile(v):
                print('Error: aenet server: file not exist: < {:} >'.format(v))
                sys.exit(-1)


    import ase
    from ase.calculators.calculator import Calculator
    from ase.neighborlist import NeighborList
    from core import ANNPotentials
    import multiprocessing


    if DEBUG:
        info = ''
        for i,m in enumerate(AENETMDS):
            info += '  => id: {:}\n'.format(i+1)
            for k,v in m.items():
                info += '     {:}: {:}\n'.format(k,v)
            info += '\n'
        with open('debuginfo.txt','a') as f:
            f.write('*'*50+'\n')
            f.write('Note: the following model files will be loaded..\n')
            f.write(info)


    def qmtxt2atoms(qmfile):
        """process qmtxt file

        Return:
            log: dict
            atoms: List[[type,x,y,z]]
        """
        log = {'nice':True, 'info':''}
        atoms = []
        with open(qmfile,'rt') as f:
            for line in f:
                ltmp = line.split()
                if len(ltmp) == 0: continue

                if len(ltmp) == 4:
                    ls = []
                    if ltmp[0] in ATOMTYPES:
                        ls.append(ltmp[0])
                    elif ltmp[0] in PT2:
                        ls.append(PT2[ltmp[0]])
                    elif ltmp[0] in PERODIC_TABLE:
                        ls.append(PERODIC_TABLE[ltmp[0]])
                    else:
                        log['nice'] = False
                        log['info'] = 'Error line: {:}'.format(line)
                        log['info'] += '  => atomtype not in PERODIC_TABLE'
                else:
                    log['nice'] = False
                    log['info'] = 'Error line: {:}'.format(line)
                    log['info'] += '  => number of indices have to be 4'

                if log['nice']:
                    try:
                        for i in ltmp[1:]:
                            ls.append(float(i))
                    except ValueError:
                        log['nice'] = False
                        log['info'] = 'Error line: {:}'.format(line)

                if log['nice']:
                    atoms.append(ls)
                else:
                    break

        if log['nice']:
            return log,atoms
        return log, ''


    class ANNCalculatorMP(Calculator):
        """Pyaenet Predictor

        Args:
            results (multiprocessing.Manager().dict()):
            atoms (List[[atomtype, x, y, z]]): will be converted to: ase.Atoms
            potfiles (dict{atomtype:file}) : training potential files
        """
        def __init__(self, results, atoms, potfiles, both=None, **kwargs):
            super().__init__(**kwargs)
            self._unit = potfiles.pop('unit',None)
            self.ann = ANNPotentials(potfiles)
            self.cutoff = self.ann.Rc_max
            self.cutoff2 = self.cutoff*self.cutoff
            self.results = results
            self.atoms = ase.Atoms(
                [i[0] for i in atoms],
                positions=[i[1:] for i in atoms]
            )
            self._cell = self.atoms.get_cell()
            # update neighbor list
            cutoffs = self.cutoff*numpy.ones(len(self.atoms))
            self.neighbors = NeighborList(cutoffs,self_interaction=False,bothways=True)
            self.neighbors.update(self.atoms)
            if both:
                self.calculate_energy_and_forces()
            else:
                self.calculate_energy()

        def calculate_energy(self):
            energy = 0.0
            atom_types = self.atoms.get_chemical_symbols()
            for i in range(len(self.atoms)):
                indices, offsets = self.neighbors.get_neighbors(i)
                itype = atom_types[i]
                icor = self.atoms.positions[i]
                jcor = numpy.empty((len(indices),3), dtype=numpy.double)
                jtype = []
                inb = 0
                for j, offset in zip(indices, offsets):
                    coo = (self.atoms.positions[j] + numpy.dot(offset,self._cell))
                    d2 = numpy.sum((coo - icor)**2)
                    if d2 <= self.cutoff2:
                        jcor[inb] = coo
                        jtype.append(atom_types[j])
                        inb += 1
                energy += self.ann.atomic_energy(icor, itype, jcor[:inb], jtype)
            if self._unit and self._unit.lower() == 'hartree': energy *= 627.509469
            self.results['energy'] = energy

        def calculate_energy_and_forces(self):
            energy = 0.0
            atom_types = self.atoms.get_chemical_symbols()
            forces = numpy.zeros((len(self.atoms), 3), dtype=numpy.double)
            for i in range(len(self.atoms)):
                indices, offsets = self.neighbors.get_neighbors(i)
                itype = atom_types[i]
                idx = i + 1
                jdx = numpy.empty(indices.size, dtype=numpy.intc)
                icor = self.atoms.positions[i]
                jcor = numpy.empty((indices.size,3), dtype=numpy.double)
                jtype = []
                inb = 0
                for j, offset in zip(indices, offsets):
                    coo = (self.atoms.positions[j] + numpy.dot(offset, self._cell))
                    d2 = numpy.sum((coo - icor)**2)
                    if d2 <= self.cutoff2:
                        jcor[inb] = coo
                        jdx[inb] = j + 1
                        jtype.append(atom_types[j])
                        inb += 1
                energy += self.ann.atomic_energy_and_forces(
                    icor, itype, idx, jcor[:inb], jtype,jdx[:inb], forces
                )
            if self._unit and self._unit.lower() == 'hartree': energy *= 627.509469
            self.results['energy'] = energy
            self.results['forces'] = forces


    def predict():
        if not os.path.isfile(QMFILE):
            print('Warning: aenet server: not exist: < {:} >'.format(QMFILE))
            with open('ELEC.energy','wt') as f: f.write('FAILED\n')
            return

        log, atoms = qmtxt2atoms(QMFILE)
        if not log['nice']:
            print(log['info'])
            with open('ELEC.energy','wt') as f: f.write('FAILED\n')
            return

        enelist = []
        for f in AENETMDS:
            d = multiprocessing.Manager().dict()
            p = multiprocessing.Process(target=ANNCalculatorMP,args=(d,atoms,f))
            p.start()
            p.join()
            enelist.append(d['energy'])

        if not enelist:
            print('Warning: puzzle: how can it be?')
            ene = 'FAILED\n'
        else:
            if DEBUG:
                with open('debuginfo.txt','a') as f:
                    for i,j in enumerate(enelist):
                        f.write('{:<3} =>  {:}\n'.format(i+1,j))
                    f.write('\n')
            ene = '{:}\n'.format(stdselection(enelist,SCALE))
        with open('ELEC.energy','wt') as f: f.write(ene)


    def accept(sock):
        conn, addr = sock.accept()
        conn.setblocking(True)
        MYSEL.register(conn, 1, data=b'aenet')


    def myconnect(key):
        sock = key.fileobj
        data = key.data
        value = sock.recv(BYTES)
        # care! command line inputs will add newline symbol
        bo = True
        if value in [b'aenet', b'aenet\n']:
            predict()
            # purpose: make client on-hold
            sock.send(b'done')
            bo = False
        # force terminate
        MYSEL.unregister(sock)
        sock.close()
        return bo


    MYSEL = selectors.DefaultSelector()
    sock = socket.socket()
    sock.bind((HOST,PORT))
    sock.listen(2)
    sock.setblocking(True)
    MYSEL.register(sock, 1, data=None)
    print('Note: aenet server is ready...')
    while True:
        events = MYSEL.select(timeout=None)
        bo = False
        for key, mask in events:
            if key.data is None:
                accept(key.fileobj)
            elif key.data == b'aenet':
                bo = myconnect(key)
            if bo:
                break
        if bo:
            print('Note: aenet server terminated')
            break

    MYSEL.close()
    print('Note: done aenet server')
    sys.exit(0)


# following is for spk progrom
import torch
from torch.utils.data import Dataset, DataLoader
from ase.db import connect
from ase.io.extxyz import read_xyz
import schnetpack as spk
from schnetpack import Properties
from schnetpack.environment import SimpleEnvironmentProvider
from schnetpack.data.atoms import get_center_of_mass


MDLIST = []
for i in glob.glob('./'+SPKMD+'*'): MDLIST.append(i)
if not MDLIST:
    print('Error: not exist: model name format < {:} >'.format(SPKMD))
    sys.exit(-1)

if DEBUG:
    with open('debuginfo.txt','a') as f:
        f.write('*'*50+'\n')
        f.write('Note: the following model files will be loaded..\n')
        for i,j in enumerate(MDLIST):
            f.write('   => {:<30}  =>  {:}\n'.format(j,i+1))
        f.write('\n')


def qmtxt2db(qmfile):
    """convert qmfile to database file

    1): convert qmtxt to extxyz file
        format
            line-1   number-of-atoms
            line-2   Properties=species:S:1:pos:R:3 energy=[energy]
            line-3   atomtype    x    y    z
            line-4   ...
    2): convert extxyz to database

    Return:
        log: dict
        db: StringIO
    """
    log = {'nice':True, 'info':''}
    fout = ''
    atnum = 0
    with open(qmfile,'rt') as f:
        for line in f:
            ltmp = line.split()
            if len(ltmp) == 0: continue

            if len(ltmp) == 4:
                if ltmp[0] in ATOMTYPES:
                    pass
                elif ltmp[0] in PT2:
                    ltmp[0] = PT2[ltmp[0]]
                elif ltmp[0] in PERODIC_TABLE:
                    ltmp[0] = PERODIC_TABLE[ltmp[0]]
                else:
                    log['nice'] = False
                    log['info'] = 'Error line: {:}'.format(line)
                    log['info'] += '  => atomtype not in PERODIC_TABLE'
            else:
                log['nice'] = False
                log['info'] = 'Error line: {:}'.format(line)
                log['info'] += '  => number of indices have to be 4'

            if log['nice']:
                new = '{:2} {:>12}   {:>12}   {:>12}\n'.format(*ltmp)
                fout += new
                atnum += 1
            else:
                break
    if not log['nice']:
        return log, ''

    fout = str(atnum) + '\nProperties=species:S:1:pos:R:3 energy=0.0\n' + fout
    fxyzio = io.StringIO(fout + '\n\n')

    fdbio = io.StringIO()
    with connect(fdbio, use_lock_file=False) as conn:
        for at in read_xyz(fxyzio):
            data = {}
            for k,v in at.info.items():
                data[k] = numpy.array([v])
            conn.write(at, data=data)
    fdbio.seek(0)

    return log, fdbio


def _convert_atoms(
    atoms,
    environment_provider=SimpleEnvironmentProvider(),
    collect_triples=False,
    centering_function=None,
    output=None):
    """
        Helper function to convert ASE atoms object to SchNetPack input format.

        Args:
            atoms (ase.Atoms): Atoms object of molecule
            environment_provider (callable): Neighbor list provider.
            collect_triples (bool, optional): Set to True if angular features are needed.
            centering_function (callable or None): Function for calculating center of
                molecule (center of mass/geometry/...). Center will be subtracted from
                positions.
            output (dict): Destination for converted atoms, if not None

    Returns:
        dict of torch.Tensor: Properties including neighbor lists and masks
            reformated into SchNetPack input format.
    """
    if output is None:
        inputs = {}
    else:
        inputs = output

    # Elemental composition
    cell = numpy.array(atoms.cell.array, dtype=numpy.float32)  # get cell array

    inputs[Properties.Z] = torch.LongTensor(atoms.numbers.astype(numpy.int32))
    positions = atoms.positions.astype(numpy.float32)
    if centering_function:
        positions -= centering_function(atoms)
    inputs[Properties.R] = torch.FloatTensor(positions)
    inputs[Properties.cell] = torch.FloatTensor(cell)

    # get atom environment
    nbh_idx, offsets = environment_provider.get_environment(atoms)

    # Get neighbors and neighbor mask
    inputs[Properties.neighbors] = torch.LongTensor(nbh_idx.astype(numpy.int32))

    # Get cells
    inputs[Properties.cell] = torch.FloatTensor(cell)
    inputs[Properties.cell_offset] = torch.FloatTensor(offsets.astype(numpy.float32))

    # If requested get neighbor lists for triples
    if collect_triples:
        nbh_idx_j, nbh_idx_k, offset_idx_j, offset_idx_k = collect_atom_triples(nbh_idx)
        inputs[Properties.neighbor_pairs_j] = torch.LongTensor(nbh_idx_j.astype(numpy.int32))
        inputs[Properties.neighbor_pairs_k] = torch.LongTensor(nbh_idx_k.astype(numpy.int32))

        inputs[Properties.neighbor_offsets_j] = torch.LongTensor(
            offset_idx_j.astype(numpy.int32)
        )
        inputs[Properties.neighbor_offsets_k] = torch.LongTensor(
            offset_idx_k.astype(numpy.int32)
        )
    return inputs


class AtomsData(Dataset):
    """
    PyTorch dataset for atomistic data. The raw data is stored in the specified
    ASE database. Use together with schnetpack.data.AtomsLoader to feed data
    to your model.

    Args:
        dbpath (str): path to directory containing database.
        subset (list, optional): indices to subset. Set to None for entire database.
        available_properties (list, optional): complete set of physical properties
            that are contained in the database.
        load_only (list, optional): reduced set of properties to be loaded
        units (list, optional): definition of units for all available properties
        environment_provider (spk.environment.BaseEnvironmentProvider): define how
            neighborhood is calculated
            (default=spk.environment.SimpleEnvironmentProvider).
        collect_triples (bool, optional): Set to True if angular features are needed.
        centering_function (callable or None): Function for calculating center of
            molecule (center of mass/geometry/...). Center will be subtracted from
            positions.
    """

    ENCODING = "utf-8"

    def __init__(
        self,
        dbpath,
        subset=None,
        available_properties=None,
        load_only=None,
        units=None,
        environment_provider=SimpleEnvironmentProvider(),
        collect_triples=False,
        centering_function=get_center_of_mass,
    ):
        # DEBUG-xz
        #
        #if not dbpath.endswith(".db"):
        #    raise AtomsDataError(
        #        "Invalid dbpath! Please make sure to add the file extension '.db' to "
        #        "your dbpath."
        #    )

        self.dbpath = dbpath

        # DEBUG-xz
        ## check if database is deprecated:
        #if self._is_deprecated():
        #    self._deprecation_update()

        self.subset = subset
        self.load_only = load_only
        self.available_properties = self.get_available_properties(available_properties)
        if load_only is None:
            self.load_only = self.available_properties
        if units is None:
            units = [1.0] * len(self.available_properties)
        self.units = dict(zip(self.available_properties, units))
        self.environment_provider = environment_provider
        self.collect_triples = collect_triples
        self.centering_function = centering_function

    def get_available_properties(self, available_properties):
        """
        Get available properties from argument or database.

        Args:
            available_properties (list or None): all properties of the dataset

        Returns:
            (list): all properties of the dataset
        """
        # DEBUG-xz
        ## use the provided list
        #if not os.path.exists(self.dbpath) or len(self) == 0:
        #    if available_properties is None:
        #        raise AtomsDataError(
        #            "Please define available_properties or set "
        #            "db_path to an existing database!"
        #        )
        #    return available_properties
        # read database properties
        with connect(self.dbpath) as conn:
            atmsrw = conn.get(1)
            db_properties = list(atmsrw.data.keys())
        # check if properties match
        if available_properties is None or set(db_properties) == set(
            available_properties
        ):
            return db_properties

        raise AtomsDataError(
            "The available_properties {} do not match the "
            "properties in the database {}!".format(available_properties, db_properties)
        )

    def create_splits(self, num_train=None, num_val=None, split_file=None):
        warnings.warn(
            "create_splits is deprecated, "
            + "use schnetpack.data.train_test_split instead",
            DeprecationWarning,
        )
        return train_test_split(self, num_train, num_val, split_file)

    def create_subset(self, idx):
        """
        Returns a new dataset that only consists of provided indices.
        Args:
            idx (numpy.ndarray): subset indices

        Returns:
            schnetpack.data.AtomsData: dataset with subset of original data
        """
        idx = numpy.array(idx)
        subidx = (
            idx if self.subset is None or len(idx) == 0 else numpy.array(self.subset)[idx]
        )
        return type(self)(
            dbpath=self.dbpath,
            subset=subidx,
            load_only=self.load_only,
            environment_provider=self.environment_provider,
            collect_triples=self.collect_triples,
            centering_function=self.centering_function,
            available_properties=self.available_properties,
        )

    def __len__(self):
        if self.subset is None:
            with connect(self.dbpath) as conn:
                return conn.count()
        return len(self.subset)

    def __getitem__(self, idx):
        at, properties = self.get_properties(idx)
        properties["_idx"] = torch.LongTensor(numpy.array([idx], dtype=numpy.int32))
        return properties

    def _subset_index(self, idx):
        # get row
        if self.subset is None:
            idx = int(idx)
        else:
            idx = int(self.subset[idx])
        return idx

    def get_atoms(self, idx):
        """
        Return atoms of provided index.

        Args:
            idx (int): atoms index

        Returns:
            ase.Atoms: atoms data

        """
        idx = self._subset_index(idx)
        with connect(self.dbpath) as conn:
            row = conn.get(idx + 1)
        at = row.toatoms()
        return at

    def get_metadata(self, key=None):
        """
        Returns an entry from the metadata dictionary of the ASE db.

        Args:
            key: Name of metadata entry. Return full dict if `None`.

        Returns:
            value: Value of metadata entry or full metadata dict, if key is `None`.

        """
        with connect(self.dbpath) as conn:
            if key is None:
                return conn.metadata
            if key in conn.metadata.keys():
                return conn.metadata[key]
        return None

    def set_metadata(self, metadata=None, **kwargs):
        """
        Sets the metadata dictionary of the ASE db.

        Args:
            metadata (dict): dictionary of metadata for the ASE db
            kwargs: further key-value pairs for convenience
        """

        # merge all metadata
        if metadata is not None:
            kwargs.update(metadata)

        with connect(self.dbpath) as conn:
            conn.metadata = kwargs

    def update_metadata(self, data):
        with connect(self.dbpath) as conn:
            metadata = conn.metadata
        metadata.update(data)
        self.set_metadata(metadata)

    def _add_system(self, conn, atoms, **properties):
        data = {}

        # add available properties to database
        for pname in self.available_properties:
            try:
                data[pname] = properties[pname]
            except:
                raise AtomsDataError("Required property missing:" + pname)

        conn.write(atoms, data=data)

    def add_system(self, atoms, **properties):
        """
        Add atoms data to the dataset.

        Args:
            atoms (ase.Atoms): system composition and geometry
            **properties: properties as key-value pairs. Keys have to match the
                `available_properties` of the dataset.

        """
        with connect(self.dbpath) as conn:
            self._add_system(conn, atoms, **properties)

    def add_systems(self, atoms_list, property_list):
        """
        Add atoms data to the dataset.

        Args:
            atoms_list (list of ase.Atoms): system composition and geometry
            property_list (list): Properties as list of key-value pairs in the same
                order as corresponding list of `atoms`.
                Keys have to match the `available_properties` of the dataset.

        """
        with connect(self.dbpath) as conn:
            for at, prop in zip(atoms_list, property_list):
                self._add_system(conn, at, **prop)

    def get_properties(self, idx):
        """
        Return property dictionary at given index.

        Args:
            idx (int): data index

        Returns:

        """
        idx = self._subset_index(idx)
        with connect(self.dbpath) as conn:
            row = conn.get(idx + 1)
        at = row.toatoms()
        # extract properties
        properties = {}
        for pname in self.load_only:
            properties[pname] = torch.FloatTensor(row.data[pname].copy())

        # extract/calculate structure
        properties = _convert_atoms(
            at,
            environment_provider=self.environment_provider,
            collect_triples=self.collect_triples,
            centering_function=self.centering_function,
            output=properties,
        )
        return at, properties

    def _get_atomref(self, property):
        """
        Returns single atom reference values for specified `property`.

        Args:
            property (str): property name

        Returns:
            list: list of atomrefs
        """
        labels = self.get_metadata("atref_labels")
        if labels is None:
            return None

        col = [i for i, l in enumerate(labels) if l == property]
        assert len(col) <= 1

        if len(col) == 1:
            col = col[0]
            atomref = numpy.array(self.get_metadata("atomrefs"))[:, col : col + 1]
        else:
            atomref = None

        return atomref

    def get_atomref(self, properties):
        """
        Return multiple single atom reference values as a dictionary.

        Args:
            properties (list or str): Desired properties for which the atomrefs are
                calculated.

        Returns:
            dict: atomic references
        """
        if type(properties) is not list:
            properties = [properties]
        return {p: self._get_atomref(p) for p in properties}

    def _is_deprecated(self):
        """
        Check if database is deprecated.

        Returns:
            (bool): True if ase db is deprecated.
        """
        # check if db exists
        if not os.path.exists(self.dbpath):
            return False

        # get properties of first atom
        with connect(self.dbpath) as conn:
            data = conn.get(1).data

        # check byte style deprecation
        if True in [pname.startswith("_dtype_") for pname in data.keys()]:
            return True
        # fallback for properties stored directly in the row
        if True in [type(val) != numpy.ndarray for val in data.values()]:
            return True

        return False

    def _deprecation_update(self):
        """
        Update deprecated database to a valid ase database.
        """
        warnings.warn(
            "The database is deprecated and will be updated automatically. "
            "The old database is moved to {}.deprecated!".format(self.dbpath)
        )

        # read old database
        atoms_list, properties_list = spk.utils.read_deprecated_database(self.dbpath)
        metadata = self.get_metadata()

        # move old database
        os.rename(self.dbpath, self.dbpath + ".deprecated")

        # write updated database
        self.set_metadata(metadata=metadata)
        with connect(self.dbpath) as conn:
            for atoms, properties in tqdm(
                zip(atoms_list, properties_list),
                "Updating new database",
                total=len(atoms_list),
            ):
                conn.write(atoms, data=properties)


def _collate_aseatoms(examples):
    """
    Build batch from systems and properties & apply padding

    Args:
        examples (list):

    Returns:
        dict[str->torch.Tensor]: mini-batch of atomistic systems
    """
    properties = examples[0]
    # initialize maximum sizes
    max_size = {
        prop: numpy.array(val.size(), dtype=numpy.int32) for prop, val in properties.items()
    }
    # get maximum sizes
    for properties in examples[1:]:
        for prop, val in properties.items():
            max_size[prop] = numpy.maximum(
                max_size[prop], numpy.array(val.size(), dtype=numpy.int32)
            )
    # initialize batch
    batch = {
        p: torch.zeros(len(examples), *[int(ss) for ss in size]).type(
            examples[0][p].type()
        )
        for p, size in max_size.items()
    }
    has_atom_mask = Properties.atom_mask in batch.keys()
    has_neighbor_mask = Properties.neighbor_mask in batch.keys()

    if not has_neighbor_mask:
        batch[Properties.neighbor_mask] = torch.zeros_like(
            batch[Properties.neighbors]
        ).float()
    if not has_atom_mask:
        batch[Properties.atom_mask] = torch.zeros_like(batch[Properties.Z]).float()

    # If neighbor pairs are requested, construct mask placeholders
    # Since the structure of both idx_j and idx_k is identical
    # (not the values), only one cutoff mask has to be generated
    if Properties.neighbor_pairs_j in properties:
        batch[Properties.neighbor_pairs_mask] = torch.zeros_like(
            batch[Properties.neighbor_pairs_j]
        ).float()

    # build batch and pad
    for k, properties in enumerate(examples):
        for prop, val in properties.items():
            shape = val.size()
            s = (k,) + tuple([slice(0, d) for d in shape])
            batch[prop][s] = val

        # add mask
        if not has_neighbor_mask:
            nbh = properties[Properties.neighbors]
            shape = nbh.size()
            s = (k,) + tuple([slice(0, d) for d in shape])
            mask = nbh >= 0
            batch[Properties.neighbor_mask][s] = mask
            batch[Properties.neighbors][s] = nbh * mask.long()

        if not has_atom_mask:
            z = properties[Properties.Z]
            shape = z.size()
            s = (k,) + tuple([slice(0, d) for d in shape])
            batch[Properties.atom_mask][s] = z > 0

        # Check if neighbor pair indices are present
        # Since the structure of both idx_j and idx_k is identical
        # (not the values), only one cutoff mask has to be generated
        if Properties.neighbor_pairs_j in properties:
            nbh_idx_j = properties[Properties.neighbor_pairs_j]
            shape = nbh_idx_j.size()
            s = (k,) + tuple([slice(0, d) for d in shape])
            batch[Properties.neighbor_pairs_mask][s] = nbh_idx_j >= 0

    return batch


class AtomsLoader(DataLoader):
    r"""
    Specialized ``torch.data.DataLoader`` which uses the correct
    collate_fn for AtomsData and provides functionality for calculating mean
    and stddev.

    Arguments:
        dataset (Dataset): dataset from which to load the data.
        batch_size (int, optional): how many samples per batch to load
            (default: 1).
        shuffle (bool, optional): set to ``True`` to have the data reshuffled
            at every epoch (default: False).
        sampler (Sampler, optional): defines the strategy to draw samples from
            the dataset. If specified, ``shuffle`` must be False.
        batch_sampler (Sampler, optional): like sampler, but returns a batch of
            indices at a time. Mutually exclusive with batch_size, shuffle,
            sampler, and drop_last.
        num_workers (int, optional): how many subprocesses to use for data
            loading. 0 means that the data will be loaded in the main process.
            (default: 0)
        collate_fn (callable, optional): merges a list of samples to form a
            mini-batch (default: collate_atons).
        pin_memory (bool, optional): If ``True``, the data loader will copy
            tensors into CUDA pinned memory before returning them.
        drop_last (bool, optional): set to ``True`` to drop the last incomplete
            batch, if the dataset size is not divisible by the batch size.
            If ``False`` and the size of dataset is not divisible by the batch
            size, then the last batch will be smaller. (default: False)
        timeout (numeric, optional): if positive, the timeout value for
            collecting a batch from workers. Should always be non-negative.
            (default: 0)
        worker_init_fn (callable, optional): If not None, this will be called
            on each worker subprocess with the worker id (an int in
            ``[0, num_workers - 1]``) as input, after seeding and before data
            loading. (default: None)

    """

    def __init__(
        self,
        dataset,
        batch_size=1,
        shuffle=False,
        sampler=None,
        batch_sampler=None,
        num_workers=0,
        collate_fn=_collate_aseatoms,
        pin_memory=False,
        drop_last=False,
        timeout=0,
        worker_init_fn=None,
    ):
        super(AtomsLoader, self).__init__(
            dataset=dataset,
            batch_size=batch_size,
            shuffle=shuffle,
            sampler=sampler,
            batch_sampler=batch_sampler,
            num_workers=num_workers,
            collate_fn=collate_fn,
            pin_memory=pin_memory,
            drop_last=drop_last,
            timeout=timeout,
            worker_init_fn=worker_init_fn,
        )


DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'
print('Note: loading data on device: {:}'.format(DEVICE))
print('Note: loading models...')
LOADEDMDLIST = []
for i in MDLIST:
    md = torch.load(i,map_location=torch.device(DEVICE))
    LOADEDMDLIST.append(md)


def predict():
    if not os.path.isfile(QMFILE):
        print('Warning: pytorch server: not exist: < {:} >'.format(QMFILE))
        with open('ELEC.energy','wt') as f: f.write('FAILED\n')
        return

    log, fdbio = qmtxt2db(QMFILE)
    if not log['nice']:
        print(log['info'])
        with open('ELEC.energy','wt') as f: f.write('FAILED\n')
        return

    atoms = AtomsData(fdbio)
    loader = AtomsLoader(atoms)

    enelist = []
    for _,batch in enumerate(loader):
        batch = {k: v.to(DEVICE) for k, v in batch.items()}
        with torch.no_grad():
            for md in LOADEDMDLIST:
                pred = md(batch)
                ene = pred['energy']
                # get tensor pure value
                enelist.append(ene.item())
        # annoying part, only first entry is chosen
        break

    if not enelist:
        print('Warning: puzzle: how can it be?')
        ene = 'FAILED\n'
    else:
        if DEBUG:
            with open('debuginfo.txt','a') as f:
                for i,j in enumerate(enelist):
                    f.write('{:<3} =>  {:}\n'.format(i+1,j))
                f.write('\n')
        ene = '{:}\n'.format(stdselection(enelist,SCALE))
    with open('ELEC.energy','wt') as f: f.write(ene)


def accept(sock):
    conn, addr = sock.accept()
    conn.setblocking(True)
    MYSEL.register(conn, 1, data=b'spk')


def myconnect(key):
    sock = key.fileobj
    data = key.data
    value = sock.recv(BYTES)
    # care! command line inputs will add newline symbol
    bo = True
    if value in [b'spk', b'spk\n']:
        predict()
        # purpose: make client on-hold
        sock.send(b'done')
        bo = False
    # force terminate
    MYSEL.unregister(sock)
    sock.close()
    return bo


MYSEL = selectors.DefaultSelector()
sock = socket.socket()
sock.bind((HOST,PORT))
sock.listen(2)
sock.setblocking(True)
MYSEL.register(sock, 1, data=None)
print('Note: PyTorch server is ready...')
while True:
    events = MYSEL.select(timeout=None)
    bo = False
    for key, mask in events:
        if key.data is None:
            accept(key.fileobj)
        elif key.data == b'spk':
            bo = myconnect(key)
        if bo:
            break
    if bo:
        print('Note: PyTorch server terminated')
        break

MYSEL.close()
print('Note: done PyTorch server')
