#!/usr/bin/python
"""ncsm_vce_calc.py

To run as a script:

    $ ncsm_vce_calc.py [-F | -f[ntv]*] Aeff4 Aeff5 Aeff6 Amin
    [Amax [nhw [n1 n2 [nshell]] | n1 n2]]

In the current directory, creates a RESULTS directory in which the
Valence Cluster Expansion is performed according to the A-prescription
given by Aeff4, Aeff5, and Aeff6.

Requires Aeff4, Aeff5, Aeff6, Amin as arguments.
If 1 additional argument given,   this is assumed to be Amax.
If 2 additional arguments given, they are assumed to be Amax nhw.
If 3 additional arguments given, they are assumed to be Amax n1 n2.
If 4 additional arguments given, they are assumed to be Amax nhw n1 n2.
If 5 additional arguments given, they are assumed to be Amax nhw n1 n2 nshell.
If -F or -f precedes the arguments, force recalculation of all steps
    (NCSD, TRDENS, and VCE)
If -f[ntv]* precedes the arguments...
    If n, force recalculation of NCSD.
    If t, force recalculation of TRDENS.
    If v, force recalculation of VCE.
If no preceding argument is given, calculations that have already been done
(i.e. outfiles are present) will not be redone.
"""

from __future__ import division

import re
from os import getcwd, path, walk, mkdir, chdir, symlink, remove, link
from subprocess import call
from sys import argv

from FGetSmallerInteraction import run as truncate_interaction
from FdoVCE import run as vce_calculation
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException

# CONSTANTS
Z_NAME_MAP = {
    1: 'h_', 2: 'he', 3: 'li', 4: 'be', 5: 'b_', 6: 'c_', 7: 'n_', 8: 'o_',
    9: 'f_', 10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p_',
    16: 's_', 17: 'cl', 18: 'ar', 19: 'k_', 20: 'ca'
}
ZNAME_FMT_ALT = '%d-'
PATH_MAIN = getcwd()
PATH_TEMPLATES = path.join(PATH_MAIN, 'templates')
PATH_RESULTS = path.join(PATH_MAIN, 'results')
DNAME_FMT_NUC = '%s%d_%d_Nhw%d_%d_%d'  # name, A, Aeff
DNAME_FMT_VCE = 'vce_presc%d,%d,%d_Nhw%d_%d_%d'
# A prescription, Nhw, n1, n2
DIR_VCE = 'vce'
REGEX_TBME = 'TBME'
REGEX_EGV = 'mfdp_*\d+\.egv'
FNAME_FMT_VCE = 'A%d.int'  # A value
FNAME_FMT_NCSD_OUT = '%s%d_%d_Nhw%d_%d_%d.out'  # name, A, Aeff, Nhw, n1, n2
FNAME_FMT_TBME = 'TBMEA2srg-n3lo2.O_%d.24'  # n1
FNAME_MFDP = 'mfdp.dat'
FNAME_TRDENS_IN = 'trdens.in'
FNAME_EGV = 'mfdp.egv'
FNAME_TRDENS_OUT = 'trdens.out'
FNAME_HEFF = 'Heff_OLS.dat'
LINE_FMT_MFDP_RESTR = ' %d %-2d %d %-2d %d %-4d ! N=%d'
N_SHELL = 1
NHW = 6
N1 = 15
N2 = 15
MAX_NMAX = 15


# FUNCTIONS
def generating_a_values(n_shell):
    """Based on the given major harmonic oscillator shell, gets the 3
    A values that are used to generate the effective Hamiltonian
    :param n_shell: major oscillator shell
    """
    a0 = int((n_shell + 2) * (n_shell + 1) * n_shell / 3 * 2)
    return a0, a0 + 1, a0 + 2


def get_name(z, z_name_map=Z_NAME_MAP, alt_name=ZNAME_FMT_ALT):
    """Given the proton number, return the short element name
    :param z: number of protons
    :param z_name_map: map from proton number to abbreviated name
    :param alt_name: alternate name format if z not in z_name_map
    """
    if z in z_name_map:
        return z_name_map[z]
    else:
        return alt_name % z


def make_base_directories(a_values, presc, results_path, a_dirpaths_map):
    """Makes directories for first 3 a values if they do not exist yet
    :param a_values: Values of A for which directories are made.
    Example: If in nshell1, would make directories for He4,5,6
    :param presc: Aeff prescription for VCE expansion
    :param results_path: Path to the directory into which these base
    :param a_dirpaths_map: map from A value to directory path
    directories are put
    """
    if not path.exists(results_path):
        mkdir(results_path)
    for a, aeff in zip(a_values, presc):
        dirpath = a_dirpaths_map[a]
        if not path.exists(dirpath):
            mkdir(dirpath)


def make_mfdp_file(z, a, aeff, n_hw, n_1, n_2, path_elt, outfile_name,
                   fname_fmt_tbme=FNAME_FMT_TBME,
                   path_temp=PATH_TEMPLATES,
                   mfdp_name=FNAME_MFDP):
    """Reads the mfdp file from path_temp 
    and rewrites it into path_elt in accordance
    ith the given z, a, aeff, nhw, n1, n2, and outfile name
    :param z: Proton number
    :param a: Mass number
    :param aeff: Effective mass number for interaction
    :param n_hw: Something something something dark side
    :param n_1: Number of allowed states for single particles
    :param n_2: Number of allowed states for two particles
    :param path_elt: path to the directory into which the mfdp file is
    to be put
    :param outfile_name: name of the outfile
    :param fname_fmt_tbme: Format string for tbme filename to be formatted
    with n1
    :param path_temp: path to the template directory
    :param mfdp_name: name of the mfdp file
    """
    temp_mfdp_path = path.join(path_temp, mfdp_name)
    mfdp_path = path.join(path_elt, mfdp_name)
    replace_map = get_mfdp_replace_map(
        fname_tbme=fname_fmt_tbme % n_1,
        outfile_name=outfile_name, z=z, a=a,
        n_hw=n_hw, n_1=n_1, n_2=n_2, aeff=aeff)
    _rewrite_file(src=temp_mfdp_path, dst=mfdp_path,
                  replace_map=replace_map)


def make_mfdp_files(z, a_range, a_presc, n_hw, n_1, n_2,
                    a_dirpath_map, a_outfile_map,
                    path_temp=PATH_TEMPLATES,
                    mfdp_name=FNAME_MFDP):
    for a, aeff in zip(a_range, a_presc):
        make_mfdp_file(z=z, a=a, aeff=aeff, n_hw=n_hw, n_1=n_1, n_2=n_2,
                       path_elt=a_dirpath_map[a],
                       outfile_name=a_outfile_map[a],
                       path_temp=path_temp, mfdp_name=mfdp_name)


def get_mfdp_replace_map(fname_tbme, outfile_name, z, a, n_hw, n_1, n_2, aeff):
    n = a - z
    if a % 2 == 0:
        tot2 = 0
    else:
        tot2 = 1
    rest_lines = get_mfdp_restrictions_lines(nmax=max(n_1, n_2))
    return {'<<TBMEFILE>>': str(fname_tbme),
            '<<OUTFILE>>': str(outfile_name),
            '<<Z>>': str(z), '<<N>>': str(n),
            '<<NHW>>': str(n_hw), '<<TOT2>>': str(tot2),
            '<<N1>>': str(n_1), '<<N2>>': str(n_2),
            '<<RESTRICTIONS>>': str(rest_lines),
            '<<AEFF>>': str(aeff)}


def get_mfdp_restrictions_lines(nmax,
                                _max_allowed_nmax=MAX_NMAX,
                                _str_rest_line=LINE_FMT_MFDP_RESTR):
    lines = list()
    for n in range(min(nmax, _max_allowed_nmax) + 1):
        i = (n + 1) * (n + 2)
        lines.append(_str_rest_line % (0, i, 0, i, 0, 2 * i, n))
    return '\n'.join(lines)


def make_trdens_file(z, a, nuc_dir,
                     path_results=PATH_RESULTS,
                     path_temp=PATH_TEMPLATES,
                     trdens_name=FNAME_TRDENS_IN):
    """Reads the trdens.in file from path_temp and rewrites it 
    into path_elt in accordance with the given z, a
    :param z: proton number
    :param a: mass number
    :param nuc_dir: directory name
    :param path_results: path to the results directory
    :param path_temp: path to the templates directory
    :param trdens_name: name of the trdens file in the templates dir
    """
    src = path.join(path_temp, trdens_name)
    path_elt = path.join(path_results, nuc_dir)
    dst = path.join(path_elt, trdens_name)
    rep_map = get_trdens_replace_map(z=z, a=a)
    _rewrite_file(src=src, dst=dst, replace_map=rep_map)


def get_trdens_replace_map(z, a):
    nnn, num_states = get_num_states(z, a)
    return {'<<NNN>>': str(nnn), '<<NUMSTATES>>': str(num_states)}


def get_num_states(z, a):
    if z == 2:
        if a == 5:
            return 1, 2
        elif a == 6:
            return 2, 5
        else:
            raise UnknownNumStatesException()
    else:
        raise UnknownNumStatesException()


class UnknownNumStatesException(Exception):
    pass


def _rewrite_file(src, dst, replace_map):
    """Reads the file given by src, replaces string elements based
       on the replace map, writes the file into dst.
    """
    # read the src file
    infile = open(src, 'r')
    read_lines = infile.readlines()
    infile.close()
    # replace strings
    write_lines = list()
    for line in read_lines:
        for k, v in replace_map.iteritems():
            if k in line:
                line = line.replace(k, str(v))
        write_lines.append(line)
    # write to the dst file
    outfile = open(dst, 'w')
    outfile.writelines(write_lines)
    outfile.close()


def truncate_space(n1, n2,
                   path_elt,
                   path_temp=PATH_TEMPLATES,
                   tbme_name_regex=REGEX_TBME,
                   fname_fmt_tbme=FNAME_FMT_TBME):
    """Run the script that truncates the space by removing extraneous
    interactions from the TBME file

    :param n1: Maximum state for single particle
    :param n2: Maximum state for two particles
    :param path_elt: Path to the directory in which the resultant TBME file
    is to be put
    :param path_temp: Path to the templates directory in which the full TBME
    files resides
    :param tbme_name_regex: Regular expression that matches only the TBME file
    in the templates directory
    :param fname_fmt_tbme: Format string for the TBME file to be formatted
    with n1
    """
    w = walk(path_temp)
    dirpath, dirnames, filenames = w.next()
    for f in filenames:
        if re.match(tbme_name_regex, f) is not None:
            tbme_filename = f
            break
    else:
        raise TBMEFileNotFoundException()
    src_path = path.join(dirpath, tbme_filename)
    dst_path = path.join(path_elt, fname_fmt_tbme % n1)
    truncate_interaction(src_path, n1, n2, dst_path)


def truncate_spaces(n1, n2,
                    dirpaths, path_temp,
                    tbme_name_regex):
    """For multiple directories, perform the operation of truncate_space

    :param n1: Something
    :param n2: Something
    :param dirpaths: Paths to the destinations
    :param path_temp: Path to the source
    :param tbme_name_regex: Regex that matches the TBME file in the source
    """
    for dirpath in dirpaths:
        truncate_space(n1=n1, n2=n2, path_elt=dirpath,
                       path_temp=path_temp,
                       tbme_name_regex=tbme_name_regex)


class TBMEFileNotFoundException(Exception):
    pass


def rename_egv_file(a6_dir, egv_name_regex, next_egv_name, force):
    """Renames the egv file from its default output name to the name needed
    for running TRDENS

    :param a6_dir: Directory in which the file resides
    :param egv_name_regex: Regular expression that matches the defualt output
    name
    :param next_egv_name: Name that the file is renamed to
    :param force: If True, replaces any existing files by the name of
    next_egv_name
    """
    next_egv_path = path.join(a6_dir, next_egv_name)
    if path.lexists(next_egv_path):
        if not force:
            return 0
        else:
            remove(next_egv_path)
    dirpath, dirnames, filenames = walk(a6_dir).next()
    for f in filenames:
        if re.match(egv_name_regex, f) is not None:
            symlink(f, next_egv_path)
            break
    else:
        raise EgvFileNotFoundException()
    return 1


class EgvFileNotFoundException(Exception):
    pass


def do_ncsd(a_values, presc, a_dirpaths_map, a_outfile_map, force):
    """Run the NCSD calculations. For each A value, does the NCSD calculation
    in each of its corresponding directories.

    :param a_values: 3-tuple of A values which generate the effective
    Hamiltonian
    :param presc: 3-tuple of Aeff values which, along with the
    A values, are used to generate the effective Hamiltonian
    :param a_dirpaths_map: Map from A values to the path to the directory in
    which NCSD is to be performed
    :param a_outfile_map: Map from A values to the .out files produced by the
    calculation. If force is False, will not do the calculation if such files
    already exist
    :param force: If true, redoes the calculations even if the output files
    already exist.
    """
    main_dir = getcwd()
    for a, aeff in zip(a_values, presc):
        if force or not path.exists(path.join(a_dirpaths_map[a],
                                              a_outfile_map[a])):
            chdir(a_dirpaths_map[a])
            call(['NCSD'])
            chdir(main_dir)


def do_trdens(a6_dir, force, outfile):
    """Run the TRDENS calculation in a6_dir

    :param a6_dir: Directory in which to run the calulation
    :param force: If True, redoes the calculation even if output files already
    exist
    :param outfile: Name of the output file generated by the TRDENS calculation.
    (If force is False, will not run if outfile already exists)
    """
    outfile_path = path.join(a6_dir, outfile)
    if path.exists(outfile_path):
        if not force:
            return 0
        else:
            remove(outfile_path)
    main_dir = getcwd()
    chdir(a6_dir)
    call(['TRDENS'])
    chdir(main_dir)
    return 1


def do_vce(
        a_values, a_prescription, a_range,
        a_outfile_map, dirpath_aeff6, dirpath_vce,
        fname_fmt_vce_int, fname_heff,
        force,
):
    """Do the VCE expansion calculation for each Aeff value in aeff_range

    :param a_values: A values used to form the effective Hamiltonian
    :param a_prescription: Range of Aeff values to evaluate based on
    the effective Hamiltonian
    :param a_range: sequence of A values for generating interaction files
    Note: All generated interaction files will be the same (linked), the
    only difference is their file names, such that the shell_calc.py
    script interprets them as interactions for different A values.
    :param a_outfile_map: Map from A values to their respective NCSD output
    files
    :param dirpath_aeff6: Directory for the 3rd A value
    :param dirpath_vce: Path to the directory in which to put generated
    interaction files
    :param fname_fmt_vce_int: Filename template for generated interaction files
    :param fname_heff: Name of the effective Hamiltonian output file
    generated by the TRDENS calculation for the 3rd A value
    :param force: If True, force redoing the calculation even if
    output files already exist
    """
    he4_fname = a_outfile_map[a_values[0]]
    he5_fname = a_outfile_map[a_values[1]]
    he6_fname = path.join(dirpath_aeff6, fname_heff)
    a0 = a_range[0]
    fpath_fmt = path.join(dirpath_vce, fname_fmt_vce_int)
    fpath = fpath_fmt % a0
    if force or not path.exists(fpath):
        vce_calculation(a_prescription, fpath, he4_fname, he5_fname, he6_fname)
    if len(a_range) > 1:
        for a in a_range[1:]:
            next_fpath = fpath_fmt % a
            if path.exists(next_fpath) and force:
                remove(next_fpath)
                link(fpath, next_fpath)
            elif not path.exists(next_fpath):
                link(fpath, next_fpath)


def ncsd_vce_calculation(
        a_prescription, a_range,
        nshell=N_SHELL, nhw=NHW, n1=N1, n2=N2,
        force_ncsd=False,
        force_trdens=False,
        force_vce=False,
        force_all=False,
        _path_results=PATH_RESULTS,
        _path_temp=PATH_TEMPLATES,
        _dir_fmt_nuc=DNAME_FMT_NUC,
        _dir_fmt_vce=DNAME_FMT_VCE,
        _dir_vce=DIR_VCE,
        _fname_regex_tbme=REGEX_TBME,
        _fname_regex_egv=REGEX_EGV,
        _fname_fmt_ncsd_out=FNAME_FMT_NCSD_OUT,
        _fname_fmt_vce=FNAME_FMT_VCE,
        _fname_mfdp=FNAME_MFDP,
        _fname_trdens_in=FNAME_TRDENS_IN,
        _fname_trdens_out=FNAME_TRDENS_OUT,
        _fname_egv_final=FNAME_EGV,
        _fname_heff=FNAME_HEFF,
):
    """Valence cluster expansion calculations within NCSM

    :param a_prescription: 3-tuple containing the Aeff values for use in
    constructing the effective interaction Hamiltonian
    effective Hamiltonian
    :param a_range: sequence of A values for which to generate interaction
    files
    :param nshell: Major harmonic oscillator shell
    :param nhw: Something something something dark side
    :param n1: Something something something dark side
    :param n2: Something something something dark side
    :param _path_results: Path to the results directory
    :param _path_temp: Path to the templates directory
    :param _dir_fmt_nuc: Template string for the directory for the base
    interactions. The format parameters are a string and two integers,
    representing (name, A, Aeff). For example, for Nshell1 with a (6, 6, 6)
    prescription, one would create he4_6, he5_6, he6_6.
    :param _dir_fmt_vce: template string for the directory for the vce
    interaction files
    :param _dir_vce: Template string for the directory for the interactions
    calculated based on the effective Hamiltonian. Should accept a string
    representing the A_prescription tuple as a format argument.
    :param _fname_regex_tbme: Regular expression that matches the TBME file in
    the templates directory
    :param _fname_regex_egv: Regular expression that matches the .egv file
    outputted by the NCSD calculation (which needs to be renamed)
    :param _fname_fmt_ncsd_out: Template string for the .out file to be
    outputted by the NCSD calculation. (This should accept the same arguments
    as the dir_fmt_nuc, so a natural choice for this string would be appending
    '.out' to the dir_fmt_nuc).
    :param _fname_fmt_vce: Template string for the name to be given to the
    .int file generated by running the FdoVCE script.
    Must accept 7 integer arguments:
        Aeff4 Aeff5 Aeff6 Nhw n1 n2 Aeff6
    :param _fname_mfdp: Name of the mfdp file in the templates directory
    :param _fname_trdens_in: Name of the trdens input file in the templates
    directory
    :param _fname_egv_final: String with which to rename the .egv file outputted
    by the NCSD calculation in prepration for the TRDENS calculation.
    :param _fname_trdens_out: Name of the trdens output file (signifying that
    this operation has been performed).
    :param _fname_heff: Name of the effective OLS Hamiltonian file
    :param force_ncsd: If True, force redoing the NCSD calculations even if
    output files already exist
    :param force_trdens: If True, force redoing the TRDENS calculations even
    if output files already exist
    :param force_vce: If True, force redoing the VCE expansion calculations
    even if output files already exist
    :param force_all: If True, force redoing all calculations
    """
    # Get a values and directory names
    a_values = generating_a_values(nshell)
    z = a_values[0] / 2
    a_dirpaths_map = dict()
    for a, aeff in zip(a_values, a_prescription):
        a_dirpaths_map[a] = path.join(
            _path_results,
            _dir_fmt_nuc % (get_name(z), a, aeff, nhw, n1, n2))
    a_outfile_map = dict()
    for a, aeff in zip(a_values, a_prescription):
        a_outfile_map[a] = path.join(
            a_dirpaths_map[a],
            _fname_fmt_ncsd_out % (get_name(z), a, aeff, nhw, n1, n2))

    # Make directories for base files
    make_base_directories(
        a_values=a_values,
        presc=a_prescription,
        results_path=_path_results,
        a_dirpaths_map=a_dirpaths_map)

    # ncsd calculations: make mfdp files, perform truncation, and run NCSD
    make_mfdp_files(
        z=z, a_range=a_values,
        a_dirpath_map=a_dirpaths_map,
        a_outfile_map=a_outfile_map,
        a_presc=a_prescription, n_hw=nhw, n_1=n1, n_2=n2,
        mfdp_name=_fname_mfdp)
    truncate_spaces(
        n1=n1, n2=n2,
        dirpaths=a_dirpaths_map.values(),
        path_temp=_path_temp,
        tbme_name_regex=_fname_regex_tbme)
    do_ncsd(
        a_values=a_values, presc=a_prescription,
        a_dirpaths_map=a_dirpaths_map, a_outfile_map=a_outfile_map,
        force=force_ncsd or force_all)

    # for the 3rd a value, make trdens file and run TRDENS
    make_trdens_file(
        z=z, a=a_values[2],
        nuc_dir=a_dirpaths_map[a_values[2]],
        path_results=_path_results,
        path_temp=_path_temp,
        trdens_name=_fname_trdens_in)
    a6_dir = a_dirpaths_map[a_values[2]]
    rename_egv_file(
        a6_dir=a6_dir,
        egv_name_regex=_fname_regex_egv,
        next_egv_name=_fname_egv_final,
        force=force_trdens or force_all)
    do_trdens(
        a6_dir=a6_dir,
        outfile=_fname_trdens_out,
        force=force_trdens or force_all,)

    # do valence cluster expansion
    if not path.exists(_dir_vce):
        mkdir(_dir_vce)
    vce_dirpath = path.join(
        _path_results, _dir_vce,
        _dir_fmt_vce % tuple(a_prescription + (nhw, n1, n2))
    )
    if not path.exists(vce_dirpath):
        mkdir(vce_dirpath)
    do_vce(
        a_values=a_values,
        a_prescription=a_prescription,
        a_range=a_range,
        dirpath_aeff6=a6_dir,
        dirpath_vce=vce_dirpath,
        a_outfile_map=a_outfile_map,
        fname_heff=_fname_heff,
        fname_fmt_vce_int=_fname_fmt_vce,
        force=force_vce or force_all
    )

    return 1


def ncsd_vce_calculations(a_prescriptions, **kwargs):
    for presc in a_prescriptions:
        ncsd_vce_calculation(a_prescription=presc, **kwargs)


# SCRIPT
def _force_from_argv0(argv0):
    if len(argv0) == 0:
        return (False,) * 4
    elif 'f' in argv0:
        force_ncsd, force_trdens, force_vce = (False,) * 3
        if 'fn' in argv0:
            force_ncsd = True
        if 'ft' in argv0:
            force_trdens = True
        if 'fv' in argv0:
            force_vce = True
        force_all = not (force_ncsd or force_trdens or force_vce)
        return force_ncsd, force_trdens, force_vce, force_all
    elif 'F' in argv0:
        return (False,) * 3 + (True,)
    else:
        return (False,) * 4


if __name__ == "__main__":
    if '-' in argv[1]:
        f_ncsd, f_trdens, f_vce, f_all = _force_from_argv0(argv[1])
        user_args = argv[2:]
    else:
        f_ncsd, f_trdens, f_vce, f_all = (None,) * 4
        user_args = argv[1:]
    a_prescription0 = tuple([int(x) for x in user_args[0:3]])
    if len(user_args) == 4:
        a_range0 = [int(user_args[3])]
        ncsd_vce_calculation(
            a_prescription=a_prescription0, a_range=a_range0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all)
    elif len(user_args) == 5:
        a_range0 = list(range(int(user_args[3]), int(user_args[4])+1))
        ncsd_vce_calculation(
            a_prescription=a_prescription0, a_range=a_range0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all)
    elif len(user_args) == 6:
        a_range0 = list(range(int(user_args[3]), int(user_args[4])+1))
        nhw_0 = int(user_args[5])
        ncsd_vce_calculation(
            a_prescription=a_prescription0, a_range=a_range0, nhw=nhw_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all)
    elif len(user_args) == 7:
        a_range0 = list(range(int(user_args[3]), int(user_args[4])+1))
        n1_0, n2_0 = [int(x) for x in user_args[5:7]]
        ncsd_vce_calculation(
            a_prescription=a_prescription0, a_range=a_range0, n1=n1_0, n2=n2_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all)
    elif len(user_args) == 8:
        a_range0 = list(range(int(user_args[3]), int(user_args[4])+1))
        nhw_0, n1_0, n2_0 = [int(x) for x in user_args[5:8]]
        ncsd_vce_calculation(
            a_prescription=a_prescription0, a_range=a_range0,
            nhw=nhw_0, n1=n1_0, n2=n2_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all)
    elif len(user_args) == 9:
        a_range0 = list(range(int(user_args[3]), int(user_args[4])+1))
        nhw_0, n1_0, n2_0, nshell_0 = [int(x) for x in user_args[5:9]]
        ncsd_vce_calculation(
            a_prescription=a_prescription0, a_range=a_range0,
            nhw=nhw_0, n1=n1_0, n2=n2_0, nshell=nshell_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all)
    else:
        raise InvalidNumberOfArgumentsException(
            '%d' % (len(argv) - 1,) +
            ' is not a valid number of arguments for ncsm_vce_calc.py.' +
            'Please enter 4-9 arguments.')
