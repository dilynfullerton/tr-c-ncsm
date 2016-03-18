#!/usr/bin/python
"""ncsm_vce_calc.py

To run as a script:

    $ ncsm_vce_calc.py [-f[ntv]{0,3}]
    [ Aeff4 Aeff5 Aeff6 | [-m|-M] Ap_min Ap_max] Amin
    [Amax [nhw [n1 n2 [nshell]] | n1 n2]]

In the current directory, creates a RESULTS directory in which the
Valence Cluster Expansion is performed according to the A-prescription(s)
given by (Aeff4, Aeff5, and Aeff6) or the range of A prescriptions
specified by Ap_min and Ap_max if -m or -M precedes the arguments.

If -F or -f precedes the arguments:
    force recalculation of all steps (NCSD, TRDENS, and VCE)
Else if -f[ntv]* precedes the arguments:
    If n, force recalculation of NCSD.
    If t, force recalculation of TRDENS.
    If v, force recalculation of VCE.
Otherwise:
    Calculations are not done if outfiles are present
If -m or -M precedes the arguments:
    The first two arguments are used to determine a range of A prescriptions
    -m --> Prescriptions are all increasing length-3 combinations of
        integers in the range [Ap_min, Ap_max]
    -M --> Prescriptions are all increasing length-3 combinations of
        numbers with repetition in the range [Ap_min, Ap_max]
Otherwise:
    The first three arguments are used to explicitly express the A
        prescription
If 1 additional argument given,   this is assumed to be Amax.
If 2 additional arguments given, they are assumed to be Amax nhw.
If 3 additional arguments given, they are assumed to be Amax n1 n2.
If 4 additional arguments given, they are assumed to be Amax nhw n1 n2.
If 5 additional arguments given, they are assumed to be Amax nhw n1 n2 nshell.
"""

from __future__ import division

import re
from os import getcwd, path, walk, mkdir, chdir, symlink, remove, link
from subprocess import Popen, PIPE
from sys import argv, stdout
from math import floor

from FGetSmallerInteraction import run as truncate_interaction
from FdoVCE import run as vce_calculation
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException

# CONSTANTS
N_SHELL = 1
NHW = 6
N1 = 15
N2 = 15
MAX_NMAX = 15

# Directories
_DPATH_MAIN = getcwd()
_DPATH_TEMPLATES = path.join(_DPATH_MAIN, 'templates')
_DPATH_RESULTS = path.join(_DPATH_MAIN, 'results')
_DNAME_FMT_NUC = '%s%d_%d_Nhw%d_%d_%d'  # name, A, Aeff
_DNAME_FMT_VCE = 'vce_presc%d,%d,%d_Nhw%d_%d_%d'  # A presc, Nhw, n1, n2
_DNAME_FMT_VCE_SF = _DNAME_FMT_VCE + '_sf%f.2'
_DNAME_VCE = 'vce'

# Files
_RGX_TBME = 'TBME'
_RGX_EGV = 'mfdp_*\d+\.egv'
_FNAME_FMT_VCE = 'A%d.int'  # A value
_FNAME_FMT_NCSD_OUT = '%s%d_%d_Nhw%d_%d_%d.out'  # name, A, Aeff, Nhw, n1, n2
_FNAME_FMT_NCSD_OUT_SF = _FNAME_FMT_NCSD_OUT[:-4] + '_sf%f.2' + '.out'
_FNAME_FMT_TBME = 'TBMEA2srg-n3lo2.O_%d.24'  # n1
_FNAME_FMT_TBME_SF = _FNAME_FMT_TBME + '_sf%f.2'  # n1 scalefactor
_FNAME_MFDP = 'mfdp.dat'
_FNAME_TRDENS_IN = 'trdens.in'
_FNAME_EGV = 'mfdp.egv'
_FNAME_TRDENS_OUT = 'trdens.out'
_FNAME_HEFF = 'Heff_OLS.dat'
_LINE_FMT_MFDP_RESTR = ' %d %-2d %d %-2d %d %-4d ! N=%d'

# output
_FNAME_NCSD_STDOUT = '__stdout_ncsd__.txt'
_FNAME_NCSD_STDERR = '__stderr_ncsd__.txt'
_FNAME_TRDENS_STDOUT = '__stdout_trdens__.txt'
_FNAME_TRDENS_STDERR = '__stderr_trdens__.txt'
WIDTH_TERM = 79
WIDTH_PROGRESS_BAR = 48
STR_PROGRESS_BAR = 'Progress: %3d/%-3d '
_STR_PROG_NCSD = 'Doing NCSD calculations for (A, Aeff) pairs...'
_STR_PROG_VCE = 'Doing VCE calculations for Aeff prescriptions...'
_STR_PROG_NCSD_EX = 'Doing NCSD calculations for A=Aeff...'

# other
_Z_NAME_MAP = {
    1: 'h_', 2: 'he', 3: 'li', 4: 'be', 5: 'b_', 6: 'c_', 7: 'n_', 8: 'o_',
    9: 'f_', 10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p_',
    16: 's_', 17: 'cl', 18: 'ar', 19: 'k_', 20: 'ca'
}
_ZNAME_FMT_ALT = '%d-'


# FUNCTIONS
def _generating_a_values(n_shell):
    """Based on the given major harmonic oscillator shell, gets the 3
    A values that are used to generate the effective Hamiltonian
    :param n_shell: major oscillator shell
    """
    a_0 = int((n_shell + 2) * (n_shell + 1) * n_shell / 3 * 2)
    return a_0, a_0 + 1, a_0 + 2


def get_name(z, z_name_map=_Z_NAME_MAP, alt_name=_ZNAME_FMT_ALT):
    """Given the proton number, return the short element name
    :param z: number of protons
    :param z_name_map: map from proton number to abbreviated name
    :param alt_name: alternate name format if z not in z_name_map
    """
    if z in z_name_map:
        return z_name_map[z]
    else:
        return alt_name % z


def make_base_directories(a_values, presc, a_aeff_to_dpath_map,
                          _dpath_results=_DPATH_RESULTS):
    """Makes directories for first 3 a values if they do not exist yet
    :param a_values: Values of A for which directories are made.
    Example: If in nshell1, would make directories for He4,5,6
    :param presc: Aeff prescription for VCE expansion
    :param _dpath_results: Path to the directory into which these base
    :param a_aeff_to_dpath_map: map from A value to directory path
    directories are put
    """
    if not path.exists(_dpath_results):
        mkdir(_dpath_results)
    for a, aeff in zip(a_values, presc):
        dirpath = a_aeff_to_dpath_map[(a, aeff)]
        if not path.exists(dirpath):
            mkdir(dirpath)


def make_mfdp_file(z, a, aeff, n_hw, n_1, n_2, path_elt, outfile_name,
                   fname_fmt_tbme=_FNAME_FMT_TBME,
                   path_temp=_DPATH_TEMPLATES,
                   mfdp_name=_FNAME_MFDP):
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


def make_mfdp_files(z, a_values, a_presc, n_hw, n_1, n_2,
                    a_aeff_to_dpath_map, a_aeff_to_outfile_fpath_map,
                    path_temp=_DPATH_TEMPLATES,
                    _fname_mfdp=_FNAME_MFDP):
    for a, aeff in zip(a_values, a_presc):
        outfile_name = path.split(a_aeff_to_outfile_fpath_map[(a, aeff)])[1]
        make_mfdp_file(z=z, a=a, aeff=aeff, n_hw=n_hw, n_1=n_1, n_2=n_2,
                       path_elt=a_aeff_to_dpath_map[(a, aeff)],
                       outfile_name=outfile_name,
                       path_temp=path_temp, mfdp_name=_fname_mfdp)


def get_mfdp_replace_map(fname_tbme, outfile_name, z, a, n_hw, n_1, n_2, aeff):
    n = a - z
    par = a % 2
    if a % 2 == 0:
        tot2 = 0
    else:
        tot2 = 1
    rest_lines = get_mfdp_restrictions_lines(nmax=max(n_1, n_2))
    return {'<<TBMEFILE>>': str(fname_tbme),
            '<<OUTFILE>>': str(outfile_name),
            '<<Z>>': str(z), '<<N>>': str(n),
            '<<NHW>>': str(n_hw), '<<PAR>>': str(par), '<<TOT2>>': str(tot2),
            '<<N1>>': str(n_1), '<<N2>>': str(n_2),
            '<<RESTRICTIONS>>': str(rest_lines),
            '<<AEFF>>': str(aeff)}


def get_mfdp_restrictions_lines(nmax,
                                _max_allowed_nmax=MAX_NMAX,
                                _str_rest_line=_LINE_FMT_MFDP_RESTR):
    lines = list()
    for n in range(min(nmax, _max_allowed_nmax) + 1):
        i = (n + 1) * (n + 2)
        lines.append(_str_rest_line % (0, i, 0, i, 0, 2 * i, n))
    return '\n'.join(lines)


def make_trdens_file(z, a, nuc_dir,
                     _dpath_results=_DPATH_RESULTS,
                     _dpath_temp=_DPATH_TEMPLATES,
                     _fname_trdens_in=_FNAME_TRDENS_IN):
    """Reads the trdens.in file from path_temp and rewrites it 
    into path_elt in accordance with the given z, a
    :param z: proton number
    :param a: mass number
    :param nuc_dir: directory name
    :param _dpath_results: path to the results directory
    :param _dpath_temp: path to the templates directory
    :param _fname_trdens_in: name of the trdens file in the templates dir
    """
    src = path.join(_dpath_temp, _fname_trdens_in)
    path_elt = path.join(_dpath_results, nuc_dir)
    dst = path.join(path_elt, _fname_trdens_in)
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


def truncate_space(
        n1, n2, dpath_elt,
        _path_temp=_DPATH_TEMPLATES,
        _tbme_name_regex=_RGX_TBME,
        _fname_fmt_tbme=_FNAME_FMT_TBME
):
    """Run the script that truncates the space by removing extraneous
    interactions from the TBME file

    :param n1: Maximum state for single particle
    :param n2: Maximum state for two particles
    :param dpath_elt: Path to the directory in which the resultant TBME file
    is to be put
    :param _path_temp: Path to the templates directory in which the full TBME
    files resides
    :param _tbme_name_regex: Regular expression that matches only the TBME file
    in the templates directory
    :param _fname_fmt_tbme: Format string for the TBME file to be formatted
    with n1
    """
    w = walk(_path_temp)
    dirpath, dirnames, filenames = w.next()
    for f in filenames:
        if re.match(_tbme_name_regex, f) is not None:
            tbme_filename = f
            break
    else:
        raise TBMEFileNotFoundException()
    src_path = path.join(dirpath, tbme_filename)
    dst_path = path.join(dpath_elt, _fname_fmt_tbme % n1)
    if not path.exists(dst_path):
        truncate_interaction(src_path, n1, n2, dst_path)
    return dst_path


def truncate_spaces(
        n1, n2, dirpaths,
        _fname_fmt_tbme=_FNAME_FMT_TBME):
    """For multiple directories, perform the operation of truncate_space

    :param n1: max allowed one-particle state
    :param n2: max allowed two-particle state
    :param dirpaths: Paths to the destination directories
    :param _fname_fmt_tbme: unformatted string template for TBME file name,
    which takes a single integer (n1) as an argument
    """
    d0 = dirpaths[0]
    fpath0 = truncate_space(n1=n1, n2=n2, dpath_elt=d0)
    if len(dirpaths) > 1:
        fname_tbme = _fname_fmt_tbme % n1
        for d in dirpaths[1:]:
            dst_path = path.join(d, fname_tbme)
            if not path.exists(dst_path):
                link(fpath0, dst_path)


class TBMEFileNotFoundException(Exception):
    pass


def rename_egv_file(
        a6_dir, force, _fname_rgx_egv=_RGX_EGV, _fname_egv_next=_FNAME_EGV,
):
    """Renames the egv file from its default output name to the name needed
    for running TRDENS

    :param a6_dir: Directory in which the file resides
    :param force: If True, replaces any existing files by the name of
    next_egv_name
    :param _fname_rgx_egv: Regular expression that matches the defualt output
    name
    :param _fname_egv_next: Name that the file is renamed to
    """
    next_egv_path = path.join(a6_dir, _fname_egv_next)
    if path.lexists(next_egv_path):
        if not force:
            return 0
        else:
            remove(next_egv_path)
    dirpath, dirnames, filenames = walk(a6_dir).next()
    for f in filenames:
        if re.match(_fname_rgx_egv, f) is not None:
            symlink(f, next_egv_path)
            break
    else:
        raise EgvFileNotFoundException()
    return 1


class EgvFileNotFoundException(Exception):
    pass


def run_ncsd(dpath, fpath_outfile, force, verbose,
             _fname_stdout=_FNAME_NCSD_STDOUT,
             _fname_stderr=_FNAME_NCSD_STDERR):
    if force or not path.exists(path.join(fpath_outfile)):
        main_dir = getcwd()
        chdir(dpath)
        args = ['NCSD']
        if verbose:
            Popen(args=args)
        else:
            p = Popen(args=args, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            fout = open(_fname_stdout, 'w')
            fout.write(out)
            fout.close()
            ferr = open(_fname_stderr, 'w')
            ferr.write(err)
            ferr.close()
        chdir(main_dir)


def run_all_ncsd(a_values, presc,
                 a_aeff_to_dpath_map, a_aeff_to_outfile_fpath_map,
                 force, verbose):
    """Run the NCSD calculations. For each A value, does the NCSD calculation
    in each of its corresponding directories.

    :param a_values: 3-tuple of A values which generate the effective
    Hamiltonian
    :param presc: 3-tuple of Aeff values which, along with the
    A values, are used to generate the effective Hamiltonian
    :param a_aeff_to_dpath_map: Map from A values to the path to the
    directory in which NCSD is to be performed
    :param a_aeff_to_outfile_fpath_map: Map from A values to the .out
    files produced by the calculation. If force is False, will not do the
    calculation if such files already exist
    :param force: If true, redoes the calculations even if the output files
    already exist.
    :param verbose: if true, prints regular output of NCSD to stdout, otherwise
    this output is suppressed
    """
    for a, aeff in zip(a_values, presc):
        dpath = a_aeff_to_dpath_map[(a, aeff)]
        fpath_outfile = path.join(dpath,
                                  a_aeff_to_outfile_fpath_map[(a, aeff)])
        run_ncsd(dpath=dpath, fpath_outfile=fpath_outfile,
                 force=force, verbose=verbose)


def run_trdens(
        a6_dir, force, verbose,
        _fname_trdens_out=_FNAME_TRDENS_OUT,
        _fname_stdout=_FNAME_TRDENS_STDOUT,
        _fname_stderr=_FNAME_TRDENS_STDERR
):
    """Run the TRDENS calculation in a6_dir

    :param a6_dir: Directory in which to run the calulation
    :param _fname_trdens_out: Name of the output file generated by the TRDENS
    calculation. (If force is False, will not run if outfile already exists)
    :param force: If True, redoes the calculation even if output files
    already exist
    :param verbose: if true, prints regular output of TRDENS to stdout,
    otherwise suppresses output and writes it to
    _fname_stdout and _fname_stderr
    :param _fname_stdout: filename to which to write standard output of
    TRDENS if verbose is false
    :param _fname_stderr: filename to which to write standard error output
    of TRDENS if verbose is false
    """
    outfile_path = path.join(a6_dir, _fname_trdens_out)
    if path.exists(outfile_path):
        if not force:
            return 0
        else:
            remove(outfile_path)
    main_dir = getcwd()
    chdir(a6_dir)
    args = ['TRDENS']
    if verbose:
        Popen(args=args)
    else:
        p = Popen(args=args, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        fout = open(_fname_stdout, 'w')
        fout.write(out)
        fout.close()
        ferr = open(_fname_stderr, 'w')
        ferr.write(err)
        ferr.close()
    chdir(main_dir)
    return 1


def run_vce(
        a_values, a_prescription, a_range,
        a_aeff_to_outfile_fname_map, dirpath_aeff6, dirpath_vce, force,
        _fname_fmt_vce_int=_FNAME_FMT_VCE,
        _fname_heff=_FNAME_HEFF,
):
    """Do the VCE expansion calculation for each Aeff value in aeff_range

    :param a_values: A values used to form the effective Hamiltonian
    :param a_prescription: Range of Aeff values to evaluate based on
    the effective Hamiltonian
    :param a_range: sequence of A values for generating interaction files
    Note: All generated interaction files will be the same (linked), the
    only difference is their file names, such that the shell_calc.py
    script interprets them as interactions for different A values.
    :param a_aeff_to_outfile_fname_map: Map from A values to their
    respective NCSD output files
    :param dirpath_aeff6: Directory for the 3rd A value
    :param dirpath_vce: Path to the directory in which to put generated
    interaction files
    :param _fname_fmt_vce_int: Filename template for generated interaction
    files
    :param _fname_heff: Name of the effective Hamiltonian output file
    generated by the TRDENS calculation for the 3rd A value
    :param force: If True, force redoing the calculation even if
    output files already exist
    """
    he4_fname = a_aeff_to_outfile_fname_map[(a_values[0], a_prescription[0])]
    he5_fname = a_aeff_to_outfile_fname_map[(a_values[1], a_prescription[1])]
    he6_fname = path.join(dirpath_aeff6, _fname_heff)
    a_0 = a_range[0]
    fpath_fmt = path.join(dirpath_vce, _fname_fmt_vce_int)
    fpath = fpath_fmt % a_0
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


def _get_a_aeff_to_dpath_map(
        a_values, a_prescription, z, nhw, n1, n2,
        _dpath_results=_DPATH_RESULTS,
        _dir_fmt_nuc=_DNAME_FMT_NUC
):
    a_paths_map = dict()
    path_fmt = path.join(_dpath_results, _dir_fmt_nuc)
    for a, aeff in zip(a_values, a_prescription):
        a_paths_map[(a, aeff)] = path_fmt % (
            get_name(z), a, aeff, nhw + a % 2, n1, n2
        )
    return a_paths_map


def _get_a_aeff_to_outfile_fpath_map(
        a_values, a_prescription, z, nhw, n1, n2, a_aeff_to_dirpath_map,
        _fname_fmt_ncsd_out=_FNAME_FMT_NCSD_OUT
):
    a_outfile_map = dict()
    for a, aeff in zip(a_values, a_prescription):
        a_outfile_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)],
            _fname_fmt_ncsd_out % (get_name(z), a, aeff, nhw + a % 2, n1, n2))
    return a_outfile_map


def _print_progress(
        completed, total, end=False,
        bar_len=WIDTH_PROGRESS_BAR, total_width=WIDTH_TERM,
        text_fmt=STR_PROGRESS_BAR
):
    text = text_fmt % (floor(completed), total)
    p = completed / total
    bar_fill = int(floor(p * bar_len))
    progress_bar = '[' + '#'*bar_fill + ' '*(bar_len - bar_fill) + ']'
    sp_fill_len = total_width - len(text) - len(progress_bar)
    if sp_fill_len < 0:
        sp_fill_len = 0
    line = '\r' + text + ' '*sp_fill_len + progress_bar
    if end:
        line += '\n'
    stdout.write(line)
    stdout.flush()


def ncsd_single_calculation(
        z, a, aeff,
        nhw=NHW, n1=N1, n2=N1,
        force=False, verbose=False,
        _path_results=_DPATH_RESULTS,
        _dir_fmt_nuc=_DNAME_FMT_NUC,
        _fname_fmt_ncsd_out=_FNAME_FMT_NCSD_OUT,
):
    if a % 2 != nhw % 2:
        nhw += 1
    # get directory path
    dirpath_nuc = path.join(
        _path_results,
        _dir_fmt_nuc % (get_name(z=z), a, aeff, nhw, n1, n2))
    outfile_ncsd = path.join(
        dirpath_nuc,
        _fname_fmt_ncsd_out % (get_name(z=z), a, aeff, nhw, n1, n2))
    # make directory
    make_base_directories(
        a_values=[a], presc=[aeff],
        _dpath_results=_path_results,
        a_aeff_to_dpath_map={(a, aeff): dirpath_nuc}
    )

    # ncsm calculations: make mfdp file, perform truncation, and run NCSD
    make_mfdp_files(
        z=z, a_values=[a], a_presc=[aeff],
        a_aeff_to_dpath_map={(a, aeff): dirpath_nuc},
        a_aeff_to_outfile_fpath_map={(a, aeff): outfile_ncsd},
        n_hw=nhw, n_1=n1, n_2=n2,
    )
    truncate_spaces(n1=n1, n2=n2, dirpaths=[dirpath_nuc])
    run_all_ncsd(
        a_values=[a], presc=[aeff],
        a_aeff_to_dpath_map={(a, aeff): dirpath_nuc},
        a_aeff_to_outfile_fpath_map={(a, aeff): outfile_ncsd},
        force=force, verbose=verbose,
    )


def ncsd_multiple_calculations(
        z, a_values, a_presc_list,
        nhw=NHW, n1=N1, n2=N1,
        force=False, verbose=False, progress=True,
        _str_prog_ncsd=_STR_PROG_NCSD,
):
    a_aeff_set = set()
    # prepare directories
    for ap in a_presc_list:
        a_aeff_set |= set(zip(a_values, ap))
    a_list = [a_aeff[0] for a_aeff in a_aeff_set]
    aeff_list = [a_aeff[1] for a_aeff in a_aeff_set]
    a_aeff_to_dir_map = _get_a_aeff_to_dpath_map(
        a_values=a_list, a_prescription=aeff_list, z=z, nhw=nhw, n1=n1, n2=n2
    )
    a_aeff_to_outfile_map = _get_a_aeff_to_outfile_fpath_map(
        a_values=a_list, a_prescription=aeff_list, z=z, nhw=nhw, n1=n1, n2=n2,
        a_aeff_to_dirpath_map=a_aeff_to_dir_map
    )
    make_base_directories(
        a_values=a_list, presc=aeff_list,
        a_aeff_to_dpath_map=a_aeff_to_dir_map
    )
    make_mfdp_files(
        z=z, a_values=a_list, a_presc=aeff_list,
        a_aeff_to_dpath_map=a_aeff_to_dir_map,
        a_aeff_to_outfile_fpath_map=a_aeff_to_outfile_map,
        n_hw=nhw, n_1=n1, n_2=n2,
    )
    truncate_spaces(n1=n1, n2=n2, dirpaths=a_aeff_to_dir_map.values())
    # do ncsd
    jobs_total = len(a_aeff_set)
    jobs_completed = 0
    if progress:
        print _str_prog_ncsd
    for a, aeff in sorted(a_aeff_set):
        if progress:
            _print_progress(jobs_completed, jobs_total)
        run_all_ncsd(
            a_values=[a], presc=[aeff],
            a_aeff_to_dpath_map=a_aeff_to_dir_map,
            a_aeff_to_outfile_fpath_map=a_aeff_to_outfile_map,
            force=force, verbose=verbose,
        )
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)


def ncsd_exact_calculations(z, a_range, nhw=NHW, n1=N1, n2=N2,
                            force=False, verbose=False, progress=True,
                            _str_prog_ncsd_ex=_STR_PROG_NCSD_EX):
    jobs_total = len(a_range)
    jobs_completed = 0
    if progress:
        print _str_prog_ncsd_ex
    for a in a_range:
        if progress:
            _print_progress(jobs_completed, jobs_total)
        ncsd_single_calculation(z=z, a=a, aeff=a, nhw=nhw, n1=n1, n2=n2,
                                force=force, verbose=verbose)
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)


def vce_single_calculation(
        z, a_values, a_prescription, a_range,
        nhw=NHW, n1=N1, n2=N1,
        force_trdens=False, force_vce=False, verbose=False,
        _dpath_results=_DPATH_RESULTS,
        _dname_vce=_DNAME_VCE,
        _dname_fmt_vce=_DNAME_FMT_VCE,
):
    """Valence cluster expansion
    :param z: protonn number
    :param a_values: 3-tuple of A values that form the base for constructing
    Heff
    :param a_prescription: 3-tuple of Aeff values used in place of the actual
    A values in constructing Heff
    :param a_range: sequence of A values for which effective interaction files
    are to be generated
    :param nhw: major oscillator model space truncation
    :param n1: max allowed single particle state
    :param n2: max allowed two-particle state
    :param force_trdens: if true, forces redoing of the TRDENS calculation,
    even if output file(s) are present
    :param force_vce: if true, force redoing of the valence cluster expansion
    and generation of interaction files, even if interaction files are already
    present
    :param verbose: if true, prints the regular output of TRDENS to stdout,
    otherwise suppresses output
    :param _dpath_results: Path to the results directory
    :param _dname_fmt_vce: template string for the directory for the vce
    interaction files
    :param _dname_vce: Template string for the directory for the interactions
    calculated based on the effective Hamiltonian. Should accept a string
    representing the A_prescription tuple as a format argument.
    outputted by the NCSD calculation (which needs to be renamed)
    directory
    by the NCSD calculation in prepration for the TRDENS calculation.
    this operation has been performed).
    """
    a_aeff_dir_map = _get_a_aeff_to_dpath_map(
        a_values=a_values, a_prescription=a_prescription,
        z=z, nhw=nhw, n1=n1, n2=n2
    )
    a_aeff_outfile_map = _get_a_aeff_to_outfile_fpath_map(
        a_values=a_values, a_prescription=a_prescription,
        z=z, nhw=nhw, n1=n1, n2=n2, a_aeff_to_dirpath_map=a_aeff_dir_map
    )
    # for the 3rd a value, make trdens file and run TRDENS
    a_aeff6 = (a_values[2], a_prescription[2])
    make_trdens_file(
        z=z, a=a_values[2],
        nuc_dir=a_aeff_dir_map[a_aeff6],
    )
    a6_dir = a_aeff_dir_map[a_aeff6]
    rename_egv_file(
        a6_dir=a6_dir,
        force=force_trdens)
    run_trdens(
        a6_dir=a6_dir,
        force=force_trdens, verbose=verbose)

    # do valence cluster expansion
    dpath_vce0 = path.join(_dpath_results, _dname_vce)
    if not path.exists(dpath_vce0):
        mkdir(dpath_vce0)
    vce_dirpath = path.join(
        _dpath_results, _dname_vce,
        _dname_fmt_vce % tuple(tuple(a_prescription) + (nhw, n1, n2))
    )
    if not path.exists(vce_dirpath):
        mkdir(vce_dirpath)
    run_vce(
        a_values=a_values,
        a_prescription=a_prescription,
        a_range=a_range,
        dirpath_aeff6=a6_dir,
        dirpath_vce=vce_dirpath,
        a_aeff_to_outfile_fname_map=a_aeff_outfile_map,
        force=force_vce
    )


def vce_multiple_calculations(
        z, a_values, a_presc_list, a_range, nhw, n1, n2,
        force_trdens, force_vce, verbose, progress,
        _str_prog_vce=_STR_PROG_VCE
):
    jobs_total = len(a_presc_list)
    jobs_completed = 0
    if progress:
        print _str_prog_vce
    for ap in a_presc_list:
        if progress:
            _print_progress(jobs_completed, jobs_total)
        vce_single_calculation(
            z=z, a_values=a_values,
            a_prescription=ap, a_range=a_range,
            nhw=nhw, n1=n1, n2=n2,
            force_trdens=force_trdens,
            force_vce=force_vce,
            verbose=verbose
        )
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)


def ncsd_vce_calculation(
        a_prescription, a_range,
        nshell=N_SHELL, nhw=NHW, n1=N1, n2=N2,
        force_ncsd=False, force_trdens=False, force_vce=False,
        force_all=False,
        verbose=False, progress=True,
):
    """Valence cluster expansion calculations within NCSM

    :param a_prescription: 3-tuple containing the Aeff values for use in
    constructing the effective interaction Hamiltonian
    effective Hamiltonian
    :param a_range: sequence of A values for which to generate interaction
    files
    :param verbose: if true, regular output of NCSD and TRDENS are printed
    to stdout, otherwise they are suppressed and written to a file instead
    :param progress: if true (and verbose false) display a progress bar for
    each calculation
    :param nshell: Major harmonic oscillator shell
    :param nhw: major oscillator shell model space trunctation
    :param n1: max allowed single-particle state
    :param n2: max allowed two-particle staet
    :param force_ncsd: If True, force redoing the NCSD calculations even if
    output files already exist
    :param force_trdens: If True, force redoing the TRDENS calculations even
    if output files already exist
    :param force_vce: If True, force redoing the VCE expansion calculations
    even if output files already exist
    :param force_all: If True, force redoing all calculations
    """
    # Get a values and directory names
    a_values = _generating_a_values(nshell)
    z = a_values[0] / 2

    # ncsd calculations: make mfdp files, perform truncation, and run NCSD
    ncsd_multiple_calculations(
        z=z, a_values=a_values,
        a_presc_list=[a_prescription],
        nhw=nhw, n1=n1, n2=n2,
        force=force_all or force_ncsd,
        verbose=verbose, progress=progress)

    # vce calculations
    vce_single_calculation(
        z=z, a_values=a_values,
        a_prescription=a_prescription,
        a_range=a_range,
        nhw=nhw, n1=n1, n2=n2,
        force_trdens=force_trdens or force_all,
        force_vce=force_vce or force_all, verbose=verbose
    )

    return 1


def ncsd_vce_calculations(
        a_prescription, a_range,
        nshell=N_SHELL, nhw=NHW, n1=N1, n2=N2,
        force_ncsd=False, force_trdens=False, force_vce=False,
        force_all=False,
        verbose=False, progress=True,
):
    a_values = _generating_a_values(n_shell=nshell)
    z = int(a_values[0] / 2)
    a_presc_list = list(a_prescription)
    ncsd_multiple_calculations(
        z=z, a_values=a_values,
        a_presc_list=a_presc_list,
        nhw=nhw, n1=n1, n2=n2,
        force=force_all or force_ncsd,
        verbose=verbose, progress=progress,
    )
    vce_multiple_calculations(
        z=z, a_values=a_values,
        a_presc_list=a_presc_list,
        a_range=a_range,
        nhw=nhw, n1=n1, n2=n2,
        force_trdens=force_trdens or force_all,
        force_vce=force_vce or force_all,
        verbose=verbose, progress=progress,
    )


def _combinations(sequence, r):
    if r == 0:
        yield []
    else:
        n = len(sequence)
        for i in range(n):
            m = _combinations(sequence[i+1:], r-1)
            for mi in m:
                yield [sequence[i]] + mi


def _multicombinations(sequence, r):
    if r == 0:
        yield []
    else:
        n = len(sequence)
        for i in range(n):
            m = _multicombinations(sequence[i:], r-1)
            for mi in m:
                yield [sequence[i]] + mi


# SCRIPT
def _force_from_argv0(argv0):
    force_ncsd, force_trdens, force_vce, force_all = (False,) * 4
    if 'n' in argv0:
        force_ncsd = True
    if 't' in argv0:
        force_trdens = True
    if 'v' in argv0:
        force_vce = True
    force_all = not (force_ncsd or force_trdens or force_vce)
    return force_ncsd, force_trdens, force_vce, force_all


# todo make this script handling better
if __name__ == "__main__":
    user_args = argv[1:]
    f_ncsd, f_trdens, f_vce, f_all = (False,) * 4
    multicom, com = False, False
    verbose0, progress0 = False, True
    while True:
        a0 = user_args[0]
        if re.match('^-f[ntv]{0,3}$', a0.lower()):
            f_ncsd, f_trdens, f_vce, f_all = _force_from_argv0(a0)
        elif re.match('^-m$', a0):
            com = True
        elif re.match('^-M$', a0):
            multicom = True
        elif re.match('^-v$', a0):
            verbose0, progress0 = True, False
        else:
            break
        user_args = user_args[1:]
    if com or multicom:
        ap_min, ap_max = [int(x) for x in user_args[0:2]]
        fn = ncsd_vce_calculations
        if com:
            a_prescription0 = _combinations(range(ap_min, ap_max+1), 3)
        else:
            a_prescription0 = _multicombinations(range(ap_min, ap_max+1), 3)
        other_args = user_args[2:]
    else:
        fn = ncsd_vce_calculation
        a_prescription0 = tuple([int(x) for x in user_args[0:3]])
        other_args = user_args[3:]
    if len(other_args) == 1:
        a_range0 = [int(other_args[0])]
        fn(
            a_prescription=a_prescription0, a_range=a_range0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0
        )
    elif len(other_args) == 2:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        fn(
            a_prescription=a_prescription0, a_range=a_range0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
        )
    elif len(other_args) == 3:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nhw_0 = int(other_args[2])
        fn(
            a_prescription=a_prescription0, a_range=a_range0, nhw=nhw_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
        )
    elif len(other_args) == 4:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        n1_0, n2_0 = [int(x) for x in other_args[2:4]]
        fn(
            a_prescription=a_prescription0, a_range=a_range0, n1=n1_0, n2=n2_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
        )
    elif len(other_args) == 5:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nhw_0, n1_0, n2_0 = [int(x) for x in other_args[2:5]]
        fn(
            a_prescription=a_prescription0, a_range=a_range0,
            nhw=nhw_0, n1=n1_0, n2=n2_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
        )
    elif len(other_args) == 6:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nhw_0, n1_0, n2_0, nshell_0 = [int(x) for x in other_args[2:6]]
        fn(
            a_prescription=a_prescription0, a_range=a_range0,
            nhw=nhw_0, n1=n1_0, n2=n2_0, nshell=nshell_0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
        )
    else:
        raise InvalidNumberOfArgumentsException(
            '%d' % (len(argv) - 1,) +
            ' is not a valid number of arguments for ncsm_vce_calc.py.' +
            'Please enter 3-9 arguments.')
