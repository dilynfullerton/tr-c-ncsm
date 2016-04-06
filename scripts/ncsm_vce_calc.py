#!/usr/bin/python
"""ncsm_vce_calc.py

To run as a script:

    $ ncsm_vce_calc.py [-f[ntv]{0,3}] [-v] [-s scalefactor] [-t walltime]
    [ Aeff4 Aeff5 Aeff6 | [-m|-M|-e] Ap_min Ap_max] Amin
    [Amax [nmax [n1 n2 [nshell [ncomponent]]] | n1 n2]]

In the current directory, creates a RESULTS directory in which the
Valence Cluster Expansion is performed according to the A-prescription(s)
given by (Aeff4, Aeff5, and Aeff6) or the range of A prescriptions
specified by Ap_min and Ap_max if -m or -M precedes the arguments.

-F or -f
    force recalculation of all steps (NCSD, TRDENS, and VCE)
-f[ntv]*
    n
        force recalculation of NCSD.
    t
        force recalculation of TRDENS.
    v
        force recalculation of VCE.
-m or -M or -e
    The first two arguments are used to determine a range of A prescriptions
    -m
        Prescriptions are all increasing length-3 combinations of
          integers in the range [Ap_min, Ap_max]
    -M
        Prescriptions are all increasing length-3 combinations of
          numbers with repetition in the range [Ap_min, Ap_max]
    -e
         Prescriptions are all (A, A, A) for A in the range [Ap_min, Ap_max]
-s scalefactor
    Off-diagonal valence space coupling terms in the interaction file are
      scaled by scalefactor
-t walltime
    NCSD jobs (submitted to cluster) are allotted the given amount of walltime,
      a string in the format hh:mm:ss

If 1 additional argument given,   this is Amax.
If 2 additional arguments given, they are Amax nmax.
If 3 additional arguments given, they are Amax      n1 n2.
If 4 additional arguments given, they are Amax nmax n1 n2.
If 5 additional arguments given, they are Amax nmax n1 n2 nshell.
If 6 additional arguments given, they are Amax nmax n1 n2 nshell ncomponent.
"""

from __future__ import division

import re
from os import getcwd, path, walk, mkdir, symlink, remove, link
from subprocess import Popen, PIPE
from sys import argv, stdout
from math import floor
from threading import Thread
from Queue import Queue

from FGetSmallerInteraction import run as truncate_interaction
from F_scale_interaction import scale_off_diag_outside_valence as scale_int
from FdoVCE import run as vce_calculation
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException

# CONSTANTS
N_SHELL = 1
N_COMPONENT = 2
NMAX = 6
N1 = 15
N2 = 15
MAX_NMAX = 15

# Directories
DPATH_MAIN = getcwd()
DPATH_TEMPLATES = path.join(DPATH_MAIN, 'templates')
DPATH_RESULTS = path.join(DPATH_MAIN, 'results')
DNAME_FMT_NUC = '%s%d_%d_Nhw%d_%d_%d'  # name, A, Aeff, nhw, n1, n2
DNAME_FMT_NUC_SF = DNAME_FMT_NUC + '_scale%.2f'  # scale factor
DNAME_FMT_VCE = 'vce_presc%d,%d,%d_Nmax%d_%d_%d_shell%d_dim%d'
#     A presc, Nmax, n1, n2, nshell, ncomponent
DNAME_FMT_VCE_SF = DNAME_FMT_VCE + '_scale%.2f'  # scale factor
DNAME_VCE = 'vce'

# Files
RGX_FNAME_TBME = 'TBME'
FNAME_FMT_EGV = 'mfdp_%d.egv'  # Nhw
FNAME_FMT_VCE_INT = 'A%d.int'  # A value
FNAME_FMT_NCSD_OUT = '%s%d_%d_Nhw%d_%d_%d.out'  # name, A, Aeff, Nhw, n1, n2
FNAME_FMT_NCSD_OUT_SF = FNAME_FMT_NCSD_OUT[:-4] + '_scale%.2f' + '.out'
FNAME_FMT_JOBSUB = FNAME_FMT_NCSD_OUT[:-4] + '.sh'
FNAME_FMT_JOBSUB_SF = FNAME_FMT_NCSD_OUT_SF[:-4] + '.sh'
# todo ^ implement
FNAME_FMT_TBME = 'TBMEA2srg-n3lo2.O_%d.24'  # n1
FNAME_FMT_TBME_SF = FNAME_FMT_TBME + '_sf%.2f'  # n1 scalefactor
FNAME_TMP_MFDP = 'mfdp.dat'
FNAME_TMP_TRDENS_IN = 'trdens.in'
FNAME_TMP_JOBSUB = 'job.sh'
LNAME_EGV = 'mfdp.egv'
FNAME_TRDENS_IN = 'trdens.in'
FNAME_TRDENS_OUT = 'trdens.out'
FNAME_HEFF = 'Heff_OLS.dat'
LINE_FMT_MFDP_RESTR = ' %d %-2d %d %-2d %d %-4d ! N=%d'

# output
FNAME_NCSD_STDOUT = '__stdout_ncsd__.txt'
FNAME_NCSD_STDERR = '__stderr_ncsd__.txt'
FNAME_TRDENS_STDOUT = '__stdout_trdens__.txt'
FNAME_TRDENS_STDERR = '__stderr_trdens__.txt'
FNAME_QSUB_STDOUT = '__stdout_qsub__.txt'
FNAME_QSUB_STDERR = '__stderr_qsub__.txt'
STR_PROG_NCSD = '  Doing NCSD calculations for (A, Aeff) pairs...'
STR_PROG_VCE = '  Doing VCE calculations for Aeff prescriptions...'
STR_PROG_NCSD_EX = '  Doing NCSD calculations for A=Aeff...'
WIDTH_TERM = 79
WIDTH_PROGRESS_BAR = 48
STR_PROGRESS_BAR = '    Progress: %3d/%-3d '

# other
Z_NAME_MAP = {
    1: 'h_', 2: 'he', 3: 'li', 4: 'be', 5: 'b_', 6: 'c_', 7: 'n_', 8: 'o_',
    9: 'f_', 10: 'ne', 11: 'na', 12: 'mg', 13: 'al', 14: 'si', 15: 'p_',
    16: 's_', 17: 'cl', 18: 'ar', 19: 'k_', 20: 'ca'
}
Z_NAME_FMT_ALT = '%d-'
NCSD_NUM_STATES = 15
NCSD_NUM_ITER = 200
MAX_OPEN_THREADS = 10
NCSD_CLUSTER_WALLTIME = '01:00:00'


# FUNCTIONS
def _generating_a_values(n_shell, n_component):
    """Based on the given major harmonic oscillator shell, gets the 3
    A values that are used to generate the effective Hamiltonian
    :param n_shell: major oscillator shell
    """
    a_0 = int((n_shell + 2) * (n_shell + 1) * n_shell / 3 * n_component)
    return a_0, a_0 + 1, a_0 + 2


def _get_name(z, z_name_map=Z_NAME_MAP, alt_name=Z_NAME_FMT_ALT):
    """Given the proton number, return the short element name
    :param z: number of protons
    :param z_name_map: map from proton number to abbreviated name
    :param alt_name: alternate name format if z not in z_name_map
    """
    if z in z_name_map:
        return z_name_map[z]
    else:
        return alt_name % z


def _make_base_directories(
        a_values, presc, a_aeff_to_dpath_map, dpath_results=DPATH_RESULTS):
    """Makes directories for first 3 a values if they do not exist yet
    :param a_values: Values of A for which directories are made.
    Example: If in nshell1, would make directories for He4,5,6
    :param presc: Aeff prescription for VCE expansion
    :param dpath_results: Path to the directory into which these base
    :param a_aeff_to_dpath_map: map from A value to directory path
    directories are put
    """
    if not path.exists(dpath_results):
        mkdir(dpath_results)
    for a, aeff in zip(a_values, presc):
        dirpath = a_aeff_to_dpath_map[(a, aeff)]
        if not path.exists(dirpath):
            mkdir(dirpath)


def _make_vce_directories(
        a_prescription, nmax, n1, n2, nshell, ncomponent,
        dpath_results=DPATH_RESULTS, dname_vce=DNAME_VCE, scalefactor=None
):
    dpath_vce0 = path.join(dpath_results, dname_vce)
    if not path.exists(dpath_vce0):
        mkdir(dpath_vce0)
    fmt = tuple(tuple(a_prescription) + (nmax, n1, n2, nshell, ncomponent))
    if scalefactor is not None:
        fmt += (scalefactor,)
        dname_fmt_vce = DNAME_FMT_VCE_SF
    else:
        dname_fmt_vce = DNAME_FMT_VCE
    vce_dirpath = path.join(dpath_results, dname_vce, dname_fmt_vce % fmt)
    if not path.exists(vce_dirpath):
        mkdir(vce_dirpath)
    return vce_dirpath


def _get_mfdp_restrictions_lines(nmax, max_allowed_nmax=MAX_NMAX):
    lines = list()
    for n in range(min(nmax, max_allowed_nmax) + 1):
        i = (n + 1) * (n + 2)
        lines.append(LINE_FMT_MFDP_RESTR % (0, i, 0, i, 0, 2 * i, n))
    return '\n'.join(lines)


def _get_mfdp_replace_map(
        fname_tbme, outfile_name, z, a, n_hw, n_1, n_2, aeff,
        num_states=NCSD_NUM_STATES, num_iter=NCSD_NUM_ITER
):
    n = a - z
    par = a % 2
    if a % 2 == 0:
        tot2 = 0
    else:
        tot2 = 1
    rest_lines = _get_mfdp_restrictions_lines(nmax=max(n_1, n_2))
    return {
        '<<TBMEFILE>>': str(fname_tbme),
        '<<OUTFILE>>': str(outfile_name),
        '<<Z>>': str(z), '<<N>>': str(n),
        '<<NHW>>': str(n_hw), '<<PAR>>': str(par), '<<TOT2>>': str(tot2),
        '<<N1>>': str(n_1), '<<N2>>': str(n_2),
        '<<RESTRICTIONS>>': str(rest_lines),
        '<<NUMST>>': str(num_states),
        '<<NUMITER>>': str(num_iter),
        '<<AEFF>>': str(aeff)
    }


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


def _make_mfdp_file(
        z, a, aeff, nhw, n1, n2, dpath_elt, fname_outfile,
        dpath_temp=DPATH_TEMPLATES,
        fname_tmp_mfdp=FNAME_TMP_MFDP,
):
    """Reads the mfdp file from path_temp 
    and rewrites it into path_elt in accordance
    ith the given z, a, aeff, nhw, n1, n2, and outfile name
    :param z: Proton number
    :param a: Mass number
    :param aeff: Effective mass number for interaction
    :param nhw: Something something something dark side
    :param n1: Number of allowed states for single particles
    :param n2: Number of allowed states for two particles
    :param dpath_elt: path to the directory into which the mfdp file is
    to be put
    :param fname_outfile: name of the outfile
    :param dpath_temp: path to the template directory
    :param fname_tmp_mfdp: name of the mfdp file
    """
    temp_mfdp_path = path.join(dpath_temp, fname_tmp_mfdp)
    mfdp_path = path.join(dpath_elt, fname_tmp_mfdp)
    replace_map = _get_mfdp_replace_map(
        fname_tbme=FNAME_FMT_TBME % n1,
        outfile_name=fname_outfile, z=z, a=a,
        n_hw=nhw, n_1=n1, n_2=n2, aeff=aeff
    )
    _rewrite_file(src=temp_mfdp_path, dst=mfdp_path, replace_map=replace_map)


def _make_mfdp_files(
        a_list, aeff_list, nhw_list, z, n_1, n_2,
        a_aeff_to_dpath_map, a_aeff_to_outfile_fpath_map,
        dpath_temp=DPATH_TEMPLATES, fname_tmp_mfdp=FNAME_TMP_MFDP
):
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        outfile_name = path.split(a_aeff_to_outfile_fpath_map[(a, aeff)])[1]
        _make_mfdp_file(
            z=z, a=a, aeff=aeff, nhw=nhw, n1=n_1, n2=n_2,
            dpath_elt=a_aeff_to_dpath_map[(a, aeff)], dpath_temp=dpath_temp,
            fname_outfile=outfile_name, fname_tmp_mfdp=fname_tmp_mfdp
        )


class UnknownNumStatesException(Exception):
    pass


def _get_num_states(z, a):
    if z == 2:
        if a == 4:
            return 0, 1  # todo this is a guess, is it correct?
        elif a == 5:
            return 1, 2
        elif a == 6:
            return 2, 5
        else:
            raise UnknownNumStatesException()
    else:
        raise UnknownNumStatesException()


def _get_trdens_replace_map(z, a):
    nnn, num_states = _get_num_states(z, a)
    return {'<<NNN>>': str(nnn), '<<NUMSTATES>>': str(num_states)}


def _make_trdens_file(
        z, a, nuc_dir,
        dpath_results=DPATH_RESULTS, dpath_temp=DPATH_TEMPLATES,
        fname_tmp_trdens_in=FNAME_TMP_TRDENS_IN
):
    """Reads the trdens.in file from path_temp and rewrites it 
    into path_elt in accordance with the given z, a
    :param z: proton number
    :param a: mass number
    :param nuc_dir: directory name
    :param dpath_results: path to the results directory
    :param dpath_temp: path to the templates directory
    :param fname_tmp_trdens_in: name of the trdens file in the templates dir
    """
    src = path.join(dpath_temp, fname_tmp_trdens_in)
    path_elt = path.join(dpath_results, nuc_dir)
    dst = path.join(path_elt, FNAME_TRDENS_IN)
    rep_map = _get_trdens_replace_map(z=z, a=a)
    _rewrite_file(src=src, dst=dst, replace_map=rep_map)


class TbmeFileNotFoundException(Exception):
    pass


def _truncate_space(
        nshell, n1, n2, dpath_elt, scalefactor,
        dpath_templates=DPATH_TEMPLATES,
):
    """Run the script that truncates the space by removing extraneous
    interactions from the TBME file
    :param n1: Maximum state for single particle
    :param n2: Maximum state for two particles
    :param dpath_elt: Path to the directory in which the resultant TBME file
    is to be put
    """
    w = walk(dpath_templates)
    dirpath, dirnames, filenames = w.next()
    for f in filenames:
        if re.match(RGX_FNAME_TBME, f) is not None:
            tbme_filename = f
            break
    else:
        raise TbmeFileNotFoundException()
    src_path = path.join(dirpath, tbme_filename)
    dst_path = path.join(dpath_elt, FNAME_FMT_TBME % n1)
    if not path.exists(dst_path):
        truncate_interaction(src_path, n1, n2, dst_path)
        if scalefactor is not None:
            next_dst_path = path.join(
                dpath_elt, FNAME_FMT_TBME_SF % (n1, scalefactor))
            scale_int(
                src=dst_path, dst=next_dst_path,
                nshell=nshell, scalefactor=scalefactor
            )
            remove(dst_path)
            dst_path = next_dst_path
    return dst_path


def _get_job_replace_map(walltime):
    return {'<<WALLTIME>>': str(walltime)}


def _make_job_submit_file(
        dst_fpath, walltime,
        dpath_temp=DPATH_TEMPLATES, fname_tmp_jobsub=FNAME_TMP_JOBSUB
):
    src_fpath = path.join(dpath_temp, fname_tmp_jobsub)
    rep_map = _get_job_replace_map(walltime=walltime)
    _rewrite_file(src=src_fpath, dst=dst_fpath, replace_map=rep_map)


def _make_job_submit_files(
        a_list, aeff_list, a_aeff_to_jobsub_fpath_map, walltime,
):
    a = a_list.pop()
    aeff = aeff_list.pop()
    dst = a_aeff_to_jobsub_fpath_map[(a, aeff)]
    _make_job_submit_file(dst_fpath=dst, walltime=walltime)
    if len(a_list) > 0:
        src = dst
        for a, aeff in zip(a_list, aeff_list):
            dst = a_aeff_to_jobsub_fpath_map[(a, aeff)]
            if path.exists(dst):
                remove(dst)
            link(src, dst)


def _truncate_spaces(nshell, n1, n2, dirpaths, scalefactor):
    """For multiple directories, perform the operation of truncate_space
    :param n1: max allowed one-particle state
    :param n2: max allowed two-particle state
    :param dirpaths: Paths to the destination directories
    """
    d0 = dirpaths[0]
    fpath0 = _truncate_space(
        nshell=nshell, n1=n1, n2=n2, dpath_elt=d0, scalefactor=scalefactor)
    if len(dirpaths) > 1:
        fname_tbme = FNAME_FMT_TBME % n1
        for d in dirpaths[1:]:
            dst_path = path.join(d, fname_tbme)
            if not path.exists(dst_path):
                link(fpath0, dst_path)


class EgvFileNotFoundException(Exception):
    pass


def _rename_egv_file(a6_dir, nhw, a6, force):
    """Renames the egv file from its default output name to the name needed
    for running TRDENS

    :param a6_dir: Directory in which the file resides
    :param nhw: major oscillator model space truncation
    :param a6: third A value, for which Heff is generated
    :param force: If True, replaces any existing files by the name of
    next_egv_name
    """
    if nhw % 2 != a6 % 2:
        nhw += 1
    next_egv_path = path.join(a6_dir, LNAME_EGV)
    if path.lexists(next_egv_path):
        if not force:
            return 0
        else:
            remove(next_egv_path)
    dirpath, dirnames, filenames = walk(a6_dir).next()
    fname_egv = FNAME_FMT_EGV % nhw
    for f in filenames:
        if f == fname_egv:
            symlink(f, next_egv_path)
            break
    else:
        raise EgvFileNotFoundException('File not found: %s' % fname_egv)
    return 1


def _run_ncsd(
        dpath, fpath_egv, force, verbose,
        fname_stdout=FNAME_NCSD_STDOUT, fname_stderr=FNAME_NCSD_STDERR
):
    if force or not path.exists(path.join(fpath_egv)):
        args = ['NCSD']
        if verbose:
            p = Popen(args=args, cwd=dpath)
            p.wait()
        else:
            p = Popen(args=args, cwd=dpath, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            if len(out) > 0:
                fout = open(path.join(dpath, fname_stdout), 'w')
                fout.write(out)
                fout.close()
            if len(err) > 0:
                ferr = open(path.join(dpath, fname_stderr), 'w')
                ferr.write(err)
                ferr.close()


def _run_trdens(
        a6_dir, force, verbose,
        fname_stdout=FNAME_TRDENS_STDOUT, fname_stderr=FNAME_TRDENS_STDERR
):
    """Run the TRDENS calculation in a6_dir
    :param a6_dir: Directory in which to run the calulation
    :param force: If True, redoes the calculation even if output files
    already exist
    :param verbose: if true, prints regular output of TRDENS to stdout,
    otherwise suppresses output and writes it to
    _fname_stdout and _fname_stderr
    :param fname_stdout: filename to which to write standard output of
    TRDENS if verbose is false
    :param fname_stderr: filename to which to write standard error output
    of TRDENS if verbose is false
    """
    outfile_path = path.join(a6_dir, FNAME_TRDENS_OUT)
    if path.exists(outfile_path):
        if not force:
            return 0
        else:
            remove(outfile_path)
    args = ['TRDENS']
    if verbose:
        p = Popen(args=args, cwd=a6_dir)
        p.wait()
    else:
        p = Popen(args=args, cwd=a6_dir, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if len(out) > 0:
            fout = open(path.join(a6_dir, fname_stdout), 'w')
            fout.write(out)
            fout.close()
        if len(err) > 0:
            ferr = open(path.join(a6_dir, fname_stderr), 'w')
            ferr.write(err)
            ferr.close()
    return 1


def _run_vce(
        a_values, a_prescription, a_range,
        a_aeff_to_outfile_fpath_map, dirpath_aeff6, dirpath_vce, force,
):
    """Do the VCE expansion calculation for each Aeff value in aeff_range

    :param a_values: A values used to form the effective Hamiltonian
    :param a_prescription: Range of Aeff values to evaluate based on
    the effective Hamiltonian
    :param a_range: sequence of A values for generating interaction files
    Note: All generated interaction files will be the same (linked), the
    only difference is their file names, such that the shell_calc.py
    script interprets them as interactions for different A values.
    :param a_aeff_to_outfile_fpath_map: Map from A values to their
    respective NCSD output files
    :param dirpath_aeff6: Directory for the 3rd A value
    :param dirpath_vce: Path to the directory in which to put generated
    interaction files
    :param force: If True, force redoing the calculation even if
    output files already exist
    """
    he4_fpath = a_aeff_to_outfile_fpath_map[(a_values[0], a_prescription[0])]
    he5_fpath = a_aeff_to_outfile_fpath_map[(a_values[1], a_prescription[1])]
    he6_fpath = path.join(dirpath_aeff6, FNAME_HEFF)
    if not path.exists(he4_fpath):
        raise NcsdOutfileNotFoundException(
            'NCSD outfile not found: %s' % he4_fpath)
    elif not path.exists(he5_fpath):
        raise NcsdOutfileNotFoundException(
            'NCSD outfile not found: %s' % he5_fpath)
    elif not path.exists(he6_fpath):
        raise NcsdOutfileNotFoundException(
            'NCSD outfile not found: %s' % he6_fpath)
    a_0 = a_range[0]
    fpath_fmt = path.join(dirpath_vce, FNAME_FMT_VCE_INT)
    fpath = fpath_fmt % a_0
    if force or not path.exists(fpath):
        vce_calculation(a_prescription, fpath, he4_fpath, he5_fpath, he6_fpath)
    if len(a_range) > 1:
        for a in a_range[1:]:
            next_fpath = fpath_fmt % a
            if path.exists(next_fpath) and force:
                remove(next_fpath)
                link(fpath, next_fpath)
            elif not path.exists(next_fpath):
                link(fpath, next_fpath)


def _get_a_aeff_to_dpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2,
        dpath_results=DPATH_RESULTS, scalefactor=None
):
    a_paths_map = dict()
    if scalefactor is not None:
        dname_fmt_nuc = DNAME_FMT_NUC_SF
    else:
        dname_fmt_nuc = DNAME_FMT_NUC
    path_fmt = path.join(dpath_results, dname_fmt_nuc)
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_paths_map[(a, aeff)] = path_fmt % args
    return a_paths_map


def _get_a_aeff_to_outfile_fpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2, a_aeff_to_dirpath_map,
        scalefactor=None
):
    a_aeff_outfile_map = dict()
    if scalefactor is not None:
        fname_fmt = FNAME_FMT_NCSD_OUT_SF
    else:
        fname_fmt = FNAME_FMT_NCSD_OUT
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_aeff_outfile_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], fname_fmt % args)
    return a_aeff_outfile_map


def _get_a_aeff_to_egv_fpath_map(
        a_list, aeff_list, nhw_list, a_aeff_to_dirpath_map,
):
    a_aeff_to_egv_map = dict()
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        a_aeff_to_egv_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], FNAME_FMT_EGV % nhw)
    return a_aeff_to_egv_map


def _get_a_aeff_to_jobsub_fpath_map(
        a_list, aeff_list, nhw_list, z, n1, n2, a_aeff_to_dirpath_map,
        scalefactor=None,
):
    a_aeff_jobsub_map = dict()
    if scalefactor is not None:
        fname_fmt = FNAME_FMT_JOBSUB_SF
    else:
        fname_fmt = FNAME_FMT_JOBSUB
    for a, aeff, nhw in zip(a_list, aeff_list, nhw_list):
        args = (_get_name(z), a, aeff, nhw, n1, n2)
        if scalefactor is not None:
            args += (scalefactor,)
        a_aeff_jobsub_map[(a, aeff)] = path.join(
            a_aeff_to_dirpath_map[(a, aeff)], fname_fmt % args)
    return a_aeff_jobsub_map


def _print_progress(
        completed, total, end=False,
        bar_len=WIDTH_PROGRESS_BAR, total_width=WIDTH_TERM,
        text_fmt=STR_PROGRESS_BAR
):
    if total > 0:
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


def _prepare_directories(
        a_list, aeff_list, nhw_list, z, n1, n2, nshell, scalefactor,
        cluster_submit=False, walltime=None, progress=False):
    a_aeff_to_dir_map = _get_a_aeff_to_dpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        z=z, n1=n1, n2=n2, scalefactor=scalefactor,
    )
    a_aeff_to_outfile_map = _get_a_aeff_to_outfile_fpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        z=z, n1=n1, n2=n2, a_aeff_to_dirpath_map=a_aeff_to_dir_map,
        scalefactor=scalefactor,
    )
    a_aeff_to_egvfile_map = _get_a_aeff_to_egv_fpath_map(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        a_aeff_to_dirpath_map=a_aeff_to_dir_map
    )
    if progress:
        print '  Making directories...'
    _make_base_directories(
        a_values=a_list, presc=aeff_list, a_aeff_to_dpath_map=a_aeff_to_dir_map
    )
    if progress:
        print '  Writing mfdp files...'
    _make_mfdp_files(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        a_aeff_to_dpath_map=a_aeff_to_dir_map,
        a_aeff_to_outfile_fpath_map=a_aeff_to_outfile_map,
        z=z, n_1=n1, n_2=n2,
    )
    if progress:
        print '  Truncating interaction to N1=%d N2=%d...' % (n1, n2)
    _truncate_spaces(
        nshell=nshell, n1=n1, n2=n2, dirpaths=a_aeff_to_dir_map.values(),
        scalefactor=scalefactor,
    )
    if cluster_submit:
        a_aeff_to_jobfile_map = _get_a_aeff_to_jobsub_fpath_map(
            a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
            z=z, n1=n1, n2=n2, a_aeff_to_dirpath_map=a_aeff_to_dir_map,
            scalefactor=scalefactor,
        )
        if progress:
            print '  Writing cluster submit files...'
        _make_job_submit_files(
            a_list=a_list, aeff_list=aeff_list,
            a_aeff_to_jobsub_fpath_map=a_aeff_to_jobfile_map,
            walltime=walltime,
        )
        return a_aeff_to_dir_map, a_aeff_to_egvfile_map, a_aeff_to_jobfile_map
    else:
        return a_aeff_to_dir_map, a_aeff_to_egvfile_map, dict()


def _ncsd_multiple_calculations_t(
        a_aeff_set, a_aeff_to_dpath_map, a_aeff_to_egvfile_map,
        force, progress=True, str_prog_ncsd=STR_PROG_NCSD,
        max_open_threads=MAX_OPEN_THREADS
):
    def _r(a_, aeff_):
        _run_ncsd(
            dpath=a_aeff_to_dpath_map[(a_, aeff_)],
            fpath_egv=a_aeff_to_egvfile_map[(a_, aeff_)],
            force=force, verbose=False
        )

    open_threads = Queue(maxsize=max_open_threads)
    todo_list = list(a_aeff_set)
    jobs_completed = 0
    jobs_total = len(todo_list)

    if progress:
        print str_prog_ncsd
    while len(todo_list) > 0 or not open_threads.empty():
        if progress:
            _print_progress(jobs_completed, jobs_total)
        # if room in queue, start new threads
        while len(todo_list) > 0 and not open_threads.full():
            a, aeff = todo_list.pop()
            t = Thread(target=_r, args=(a, aeff))
            open_threads.put(t)
            t.start()
        # wait for completion of first thread in queue
        t = open_threads.get()
        t.join()
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)

    if jobs_completed == jobs_total:
        return 1
    else:
        return 0


def _ncsd_multiple_calculations_s(
        a_aeff_set,
        a_aeff_to_dpath_map, a_aeff_to_egvfile_map, a_aeff_to_jobfile_map,
        force, progress,
        fname_stdout=FNAME_QSUB_STDOUT, fname_stderr=FNAME_QSUB_STDERR,
):
    submitted_jobs = 0
    if progress:
        print '  Submitting jobs...'
    for a, aeff in a_aeff_set:
        fpath_egv = a_aeff_to_egvfile_map[(a, aeff)]
        if force or not path.exists(fpath_egv):
            job = a_aeff_to_jobfile_map[(a, aeff)]
            dpath = a_aeff_to_dpath_map[(a, aeff)]
            args = ['qsub', '%s' % job]
            p = Popen(args=args, cwd=dpath, stdout=PIPE, stderr=PIPE)
            out, err = p.communicate()
            if len(out) > 0:
                fout = open(path.join(dpath, fname_stdout), 'w')
                fout.write(out)
                fout.close()
            if len(err) > 0:
                ferr = open(path.join(dpath, fname_stderr), 'w')
                ferr.write(err)
                ferr.close()
            submitted_jobs += 1
        else:
            print 'Path found: %s' % fpath_egv
    if progress:
        if submitted_jobs == 1:
            print '  1 job submitted to cluster'
        else:
            print '  %d jobs submitted to cluster' % submitted_jobs


def _ncsd_multiple_calculations(
        a_aeff_set, a_aeff_to_dpath_map, a_aeff_to_egvfile_map,
        force, verbose, progress, str_prog_ncsd=STR_PROG_NCSD
):
    # do ncsd
    jobs_total = len(a_aeff_set)
    jobs_completed = 0
    progress = progress and not verbose
    if progress:
        print str_prog_ncsd
    for a, aeff in sorted(a_aeff_set):
        if progress:
            _print_progress(jobs_completed, jobs_total)
        _run_ncsd(
            dpath=a_aeff_to_dpath_map[(a, aeff)],
            fpath_egv=a_aeff_to_egvfile_map[(a, aeff)],
            force=force, verbose=verbose
        )
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)


def ncsd_multiple_calculations(
        a_presc_list, a_values, z, nmax, a_0,
        nshell=N_SHELL,  n1=N1, n2=N1, scalefactor=None,
        force=False, verbose=False, progress=True, threading=True,
        cluster_submit=False, walltime=None,
        str_prog_ncsd=STR_PROG_NCSD,
):
    """For a given list of A prescriptions, do the NCSD calculations
    necessary for doing a valence cluster expansion
    :param a_presc_list: sequence of A prescriptions
    :param a_values: three base a values (e.g. 4, 5, 6 for p shell)
    :param z: proton number
    :param nmax: major oscillator shell model space truncation
    :param a_0: core A (e.g. 4 for p-shell, 16 for sd-shell, ...)
    :param nshell: shell (0: s, 1: p, 2: sd, ...)
    :param n1: max allowed 1-particle state
    :param n2: max allowed 2-particle state
    :param scalefactor: float value by which the off-diagonal valence
    coupling terms in the interaction are scaled; if None, no scaling is done
    :param force: if true, force doing the NCSD calculation, even if it has
    already been done
    :param verbose: if true, prints the regular NCSD output to stdout; else
    this is suppressed and written to a file instead
    :param progress: if true, show a progress bar. Note: will not show
    progress bar if verbose is true
    :param threading: if true, starts the various calculations in separate
    threads
    :param cluster_submit: if true, submits the NCSD job to cluster
    :param walltime: if cluster submit is true, this string in the format
    hh:mm:ss specifies how much time is to be allotted each NCSD
    calculation
    :param str_prog_ncsd: string to show before progress bar
    """
    a_aeff_nhw_set = set()
    # prepare directories
    for ap in a_presc_list:
        a_aeff_nhw_set |= set(zip(
            a_values, ap, [nmax + a - a_0 for a in a_values]
        ))
    a_list, aeff_list, nhw_list = list(), list(), list(),
    for a, aeff, nhw in a_aeff_nhw_set:
        a_list.append(a)
        aeff_list.append(aeff)
        nhw_list.append(nhw)
    a_aeff_maps = _prepare_directories(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        z=z, n1=n1, n2=n2, nshell=nshell, scalefactor=scalefactor,
        cluster_submit=cluster_submit, walltime=walltime, progress=progress,
    )
    a_aeff_to_dir, a_aeff_to_egv, a_aeff_to_job = a_aeff_maps
    a_aeff_set = set([(a, aeff) for a, aeff, nhw in a_aeff_nhw_set])
    if cluster_submit:
        _ncsd_multiple_calculations_s(
            a_aeff_set=a_aeff_set,
            a_aeff_to_dpath_map=a_aeff_to_dir,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            a_aeff_to_jobfile_map=a_aeff_to_job,
            progress=progress, force=force,
        )
    elif threading and len(a_aeff_nhw_set) > 1:
        _ncsd_multiple_calculations_t(
            a_aeff_set=a_aeff_set,
            a_aeff_to_dpath_map=a_aeff_to_dir,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            force=force, progress=progress,
            str_prog_ncsd=str_prog_ncsd
        )
    else:
        _ncsd_multiple_calculations(
            a_aeff_set=a_aeff_set,
            a_aeff_to_dpath_map=a_aeff_to_dir,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            force=force, verbose=verbose, progress=progress,
            str_prog_ncsd=str_prog_ncsd
        )


def ncsd_single_calculation(
        z, a, aeff, scalefactor, nhw=NMAX, n1=N1, n2=N1, nshell=N_SHELL,
        force=False, verbose=False, progress=True,
        cluster_submit=False, walltime=None,
):
    a_aeff_to_dpath, a_aeff_to_egv, a_aeff_to_job = _prepare_directories(
        a_list=[a], aeff_list=[aeff], nhw_list=[nhw],
        z=z, n1=n1, n2=n2, nshell=nshell, scalefactor=scalefactor,
        cluster_submit=cluster_submit, walltime=walltime, progress=progress,
    )[:3]
    if cluster_submit:
        return _ncsd_multiple_calculations_s(
            a_aeff_set=set([(a, aeff)]),
            a_aeff_to_dpath_map=a_aeff_to_dpath,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            a_aeff_to_jobfile_map=a_aeff_to_job,
            force=force, progress=progress,
        )
    else:
        return _ncsd_multiple_calculations(
            a_aeff_set=set([(a, aeff)]),
            a_aeff_to_dpath_map=a_aeff_to_dpath,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            force=force, progress=progress,
            verbose=verbose,
        )


def ncsd_exact_calculations(
        z, a_range, nmax=NMAX, n1=N1, n2=N2,
        nshell=N_SHELL, ncomponent=N_COMPONENT, int_scalefactor=None,
        force=False, verbose=False, progress=True,
        cluster_submit=False, walltime=None,
        str_prog_ncsd_ex=STR_PROG_NCSD_EX
):
    """For each A in a_range, does the NCSD calculation for A=Aeff
    :param z: proton number
    :param a_range: range of A values for which to do NCSD with Aeff=A
    :param nmax: major oscillator model space truncation. Note: Increased by 1
    for each successive A value
    :param n1: max allowed 1-particle state
    :param n2: max allowed 2-particle state
    :param nshell: nuclear shell (e.g. 0=s, 1=p, 2=sd, ...)
    :param ncomponent: 1=neutrons, 2=protons and neutrons
    :param int_scalefactor: float factor by which the off-diagonal valence
    coupling terms in the TBME interaction are scaled
    :param force: if true, force calculation of NCSD even if output files are
    present
    :param verbose: if true, print regular NCSD output to stdout; otherwise
    output is suppressed and written to a file instead
    :param progress: if true, show a progress bar. Note: will not be shown if
    verbose is true.
    :param cluster_submit: if true, submit job to cluster
    :param walltime: if cluster_submit is true, this string hh:mm:ss specifies
    how much wall time is to be allotted each NCSD calculation
    :param str_prog_ncsd_ex: string to show before progress bar
    """
    ncsd_multiple_calculations(
        z=z, a_values=a_range, a_presc_list=[a_range],
        nmax=nmax, n1=n1, n2=n2, nshell=nshell, scalefactor=int_scalefactor,
        a_0=_generating_a_values(n_shell=nshell, n_component=ncomponent)[0],
        cluster_submit=cluster_submit, walltime=walltime,
        force=force, verbose=verbose, progress=progress,
        str_prog_ncsd=str_prog_ncsd_ex
    )


class NcsdOutfileNotFoundException(Exception):
    pass


def vce_single_calculation(
        a_values, a_prescription, a_range, z, nmax,
        n1=N1, n2=N1, nshell=-1, ncomponent=-1, int_scalefactor=None,
        force_trdens=False, force_vce=False, verbose=False,
):
    """Valence cluster expansion
    :param z: protonn number
    :param a_values: 3-tuple of A values that form the base for constructing
    Heff
    :param a_prescription: 3-tuple of Aeff values used in place of the actual
    A values in constructing Heff
    :param a_range: sequence of A values for which effective interaction files
    are to be generated
    :param nmax: major oscillator model space truncation
    :param n1: max allowed single particle state
    :param n2: max allowed two-particle state
    :param nshell: shell number (0 = s, 1 = p, 2 = sd, ...)
    :param ncomponent: dimension (1 = neutrons, 2 = protons and neutrons)
    :param int_scalefactor: factor by which to scale the off-diagonal
    valence coupling terms of the TBME interaction
    :param force_trdens: if true, forces redoing of the TRDENS calculation,
    even if output file(s) are present
    :param force_vce: if true, force redoing of the valence cluster expansion
    and generation of interaction files, even if interaction files are already
    present
    :param verbose: if true, prints the regular output of TRDENS to stdout,
    otherwise suppresses output
    """
    a_aeff_dir_map = _get_a_aeff_to_dpath_map(
        a_list=a_values, aeff_list=a_prescription,
        nhw_list=list(range(nmax, nmax+3)), z=z, n1=n1, n2=n2
    )
    a_aeff_outfile_map = _get_a_aeff_to_outfile_fpath_map(
        a_list=a_values, aeff_list=a_prescription,
        nhw_list=list(range(nmax, nmax+3)), z=z, n1=n1, n2=n2,
        a_aeff_to_dirpath_map=a_aeff_dir_map
    )
    for f in a_aeff_outfile_map.values():
        if not path.exists(f):
            raise NcsdOutfileNotFoundException(
                'NCSD outfile not found: %s' % f)
    # for the 3rd a value, make trdens file and run TRDENS
    a_aeff6 = (a_values[2], a_prescription[2])
    _make_trdens_file(z=z, a=a_values[2], nuc_dir=a_aeff_dir_map[a_aeff6])
    a6_dirpath = a_aeff_dir_map[a_aeff6]
    try:
        _rename_egv_file(
            a6_dir=a6_dirpath, nhw=nmax+2, a6=a_values[2], force=force_trdens)
    except EgvFileNotFoundException:
        raise
    _run_trdens(a6_dir=a6_dirpath, force=force_trdens, verbose=verbose)

    # do valence cluster expansion
    vce_dirpath = _make_vce_directories(
        a_prescription=a_prescription, nmax=nmax, n1=n1, n2=n2,
        nshell=nshell, ncomponent=ncomponent, scalefactor=int_scalefactor,
    )
    try:
        _run_vce(
            a_values=a_values, a_prescription=a_prescription, a_range=a_range,
            dirpath_aeff6=a6_dirpath, dirpath_vce=vce_dirpath,
            a_aeff_to_outfile_fpath_map=a_aeff_outfile_map, force=force_vce
        )
    except NcsdOutfileNotFoundException:
        raise


def _vce_multiple_calculations_t(
        z, a_values, a_presc_list, a_range, nmax, n1, n2, nshell, ncomponent,
        force_trdens, force_vce, verbose, progress,
        int_scalefactor=None, max_open_threads=MAX_OPEN_THREADS,
        str_prog_vce=STR_PROG_VCE,
):
    error_messages = Queue()

    def _r(ap0):
        try:
            vce_single_calculation(
                z=z, a_values=a_values,
                a_prescription=ap0, a_range=a_range,
                nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
                int_scalefactor=int_scalefactor,
                force_trdens=force_trdens, force_vce=force_vce, verbose=verbose,
            )
        except EgvFileNotFoundException, e:
            error_messages.put(str(e))
        except NcsdOutfileNotFoundException, e:
            error_messages.put(str(e))

    open_threads = Queue(maxsize=max_open_threads)
    todo_list = list(a_presc_list)
    jobs_completed = 0
    jobs_total = len(todo_list)

    if progress:
        print str_prog_vce
    while len(todo_list) > 0 or not open_threads.empty():
        if progress:
            _print_progress(jobs_completed, jobs_total)
        # if room in queue, start new threads
        while len(todo_list) > 0 and not open_threads.full():
            ap = todo_list.pop()
            t = Thread(target=_r, args=(ap,))
            open_threads.put(t)
            t.start()
        # wait for completion of first thread in queue
        t = open_threads.get()
        t.join()
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)

    while not error_messages.empty():
        em = error_messages.get()
        print em

    if jobs_completed == jobs_total:
        return 1
    else:
        return 0


def _vce_multiple_calculations(
        z, a_values, a_presc_list, a_range, nmax, n1, n2, nshell, ncomponent,
        force_trdens, force_vce, verbose, progress,
        int_scalefactor=None,
        str_prog_vce=STR_PROG_VCE,
):
    jobs_total = len(a_presc_list)
    jobs_completed = 0
    progress = progress and not verbose
    if progress:
        print str_prog_vce
    error_messages = list()
    for ap in a_presc_list:
        if progress:
            _print_progress(jobs_completed, jobs_total)
        try:
            vce_single_calculation(
                z=z, a_values=a_values,
                a_prescription=ap, a_range=a_range,
                nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
                int_scalefactor=int_scalefactor,
                force_trdens=force_trdens, force_vce=force_vce, verbose=verbose,
            )
            jobs_completed += 1
        except NcsdOutfileNotFoundException, e:
            error_messages.append(str(e))
            continue
        except EgvFileNotFoundException, e:
            error_messages.append(str(e))
            continue
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)
    for em in error_messages:
        print em
    if jobs_completed == jobs_total:
        return 1
    else:
        return 0


def vce_multiple_calculations(
        a_values, a_presc_list, a_range, z, nmax, n1, n2, nshell, ncomponent,
        force_trdens, force_vce, verbose, progress, threading,
        int_scalefactor=None,
):
    """For every A prescription in a_presc_list, performs the valence
    cluster expansion to generate an interaction file. Assumes NCSD
    calculations have already been done sufficiently.
    :param z: proton number
    :param a_values: three base A values (e.g. if in p shell these are 4, 5, 6)
    :param a_presc_list: sequence of 3-tuples representing A prescriptions
    with which to do the valence cluster expansion
    :param a_range: sequence of values for which to generate interaction
    files
    :param nmax: major oscillator model space truncation
    :param n1: max allowed one-particle state
    :param n2: max allowed two-particle state
    :param nshell: shell number (0=s, 1=p, 2=sd, ...)
    :param ncomponent: number of components
    (1=neutrons, 2=protons and neutrons)
    :param int_scalefactor: factor by which TBME interaction file was scaled
    (for directory naming purposes); if None, directory will be named as
    usual
    :param force_trdens: if true, force calculation of TRDENS even if output
    file is already present
    :param force_vce: if true, force calculation of valence cluster expansion
    even if interaction output file is already present
    :param verbose: if true, regular output of TRDENS is printed to stdout;
    otherwise output is suppressed and written to a file instead
    :param progress: if true, display a progress bar. Note: if verbose is
    true, progress bar will not be displayed.
    :param threading: if true, calculations will be multi-threaded
    """
    if threading and len(a_presc_list) > 1:
        _vce_multiple_calculations_t(
            a_values=a_values, a_presc_list=a_presc_list, a_range=a_range,
            z=z, nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
            int_scalefactor=int_scalefactor,
            force_trdens=force_trdens, force_vce=force_vce,
            verbose=verbose, progress=progress,
        )
    else:
        return _vce_multiple_calculations(
            a_values=a_values, a_presc_list=a_presc_list, a_range=a_range,
            z=z, nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
            int_scalefactor=int_scalefactor,
            force_trdens=force_trdens, force_vce=force_vce,
            verbose=verbose, progress=progress,
        )


def ncsd_vce_calculations(
        a_prescriptions, a_range,
        nmax=NMAX, n1=N1, n2=N2, nshell=N_SHELL, ncomponent=N_COMPONENT,
        int_scalefactor=None,
        force_ncsd=False, force_trdens=False, force_vce=False, force_all=False,
        verbose=False, progress=True, threading=True,
        cluster_submit=False, walltime=None,
):
    """Given a sequence or generator of A prescriptions, does the NCSD/VCE
    calculation for each prescription
    :param a_prescriptions: sequence or generator of A prescription tuples
    :param a_range: sequence of A values for which to generate interaction
    files
    :param nmax: model space max oscillator shell
    :param n1: max allowed one-particle state
    :param n2: max allowed two-particle state
    :param nshell: major oscillator shell (0 = s, 1 = p, 2 = sd, ...)
    :param ncomponent: num components (1=neutrons, 2=protons and neutrons)
    :param int_scalefactor: float factor by which the off-diagonal valence
    coupling terms in the interaction are scaled
    :param force_ncsd: if true, forces recalculation of NCSD, even if output
    files already exist
    :param force_trdens: if true, forces recalculation of TRDENS, even if
    output file already exists
    :param force_vce: if true, forces recalculation of valence cluster
    expansion, even if interaction file already exists
    :param force_all: if true, forces recalculation of everything
    :param verbose: if true, prints standard output for subprocesses to
    stdout; otherwise this is suppressed and saved in a text file
    :param progress: if true, shows a progress bar. Note if verbose is true,
    progress bar will not be shown.
    :param threading: if true, calculation are multi-threaded
    :param cluster_submit: if true, NCSD calculation jobs will be submitted
    to the OpenMP cluster. NOTE: This may result in calculations remaining
    undone.
    :param walltime: if cluster_submit, this string hh:mm:ss specifies how
    much wall time is to be allotted each NCSD calculation
    """
    a_values = _generating_a_values(n_shell=nshell, n_component=ncomponent)
    z = int(a_values[0] / ncomponent)
    a_presc_list = list(a_prescriptions)
    ncsd_multiple_calculations(
        z=z, a_values=a_values, a_presc_list=a_presc_list,
        nmax=nmax, a_0=a_values[0], n1=n1, n2=n2, nshell=nshell,
        scalefactor=int_scalefactor,
        force=force_all or force_ncsd,
        verbose=verbose, progress=progress, threading=threading,
        cluster_submit=cluster_submit, walltime=walltime,
    )
    vce_multiple_calculations(
        z=z, a_values=a_values, a_presc_list=a_presc_list, a_range=a_range,
        nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
        int_scalefactor=int_scalefactor,
        force_trdens=force_trdens or force_all,
        force_vce=force_vce or force_all,
        verbose=verbose, progress=progress, threading=threading,
    )


def _exact(sequence, r):
    for s in sequence:
        yield [s] * r


def _combinations(sequence, r):
    """Given a sequence of unique elements, yields all unique length-r
    combinations of those elements WITHOUT repeats
    :param sequence: sequence of unique elements
    :param r: length of the combinations to return
    """
    if r == 0:
        yield []
    else:
        n = len(sequence)
        for i in range(n):
            m = _combinations(sequence[i+1:], r-1)
            for mi in m:
                yield [sequence[i]] + mi


def _multicombinations(sequence, r):
    """Given a sequence, yields all possible length-r multi-combinations
    of items of the sequence, where a multi-combination is
    a unique combination of elements where repeats ARE allowd
    :param sequence: sequence of unique elements from which to generate
    multi-combinations
    :param r: length of the multi-combinations to generate
    """
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
    multicom, com, exact = (False,) * 3
    verbose0, progress0 = False, True
    cluster_submit0 = False
    walltime0 = NCSD_CLUSTER_WALLTIME
    scalefactor0 = None
    while True:
        a0 = user_args[0]
        if re.match('^-f[ntv]{0,3}$', a0.lower()):
            f_ncsd, f_trdens, f_vce, f_all = _force_from_argv0(a0)
        elif a0 == '-m' or a0 == '--combinations':
            com = True
        elif a0 == '-M' or a0 == '--multicombinations':
            multicom = True
        elif a0 == '-e' or a0 == '--exact':
            exact = True
        elif a0 == '-v' or a0 == '--verbose':
            verbose0, progress0 = True, False
        elif a0 == '-s' or a0 == '--scale-int':
            user_args = user_args[1:]
            scalefactor0 = round(float(user_args[0]), 2)
        elif a0 == '-t' or a0 == '--walltime':
            user_args = user_args[1:]
            cluster_submit0 = True
            walltime0 = user_args[0]
        else:
            break
        user_args = user_args[1:]
    if com or multicom or exact:
        ap_min, ap_max = [int(x) for x in user_args[0:2]]
        ap_range = range(ap_min, ap_max+1)
        if com:
            a_prescriptions0 = _combinations(ap_range, 3)
        elif multicom:
            a_prescriptions0 = _multicombinations(ap_range, 3)
        else:
            a_prescriptions0 = _exact(ap_range, 3)
        other_args = user_args[2:]
    else:
        a_prescriptions0 = [tuple([int(x) for x in user_args[0:3]])]
        other_args = user_args[3:]
    if len(other_args) == 1:
        a_range0 = [int(other_args[0])]
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0,
            int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    elif len(other_args) == 2:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0,
            int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    elif len(other_args) == 3:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nmax_0 = int(other_args[2])
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0, nmax=nmax_0,
            int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    elif len(other_args) == 4:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        n1_0, n2_0 = [int(x) for x in other_args[2:]]
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0,
            n1=n1_0, n2=n2_0, int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    elif len(other_args) == 5:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nmax_0, n1_0, n2_0 = [int(x) for x in other_args[2:]]
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0,
            nmax=nmax_0, n1=n1_0, n2=n2_0, int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    elif len(other_args) == 6:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nmax_0, n1_0, n2_0, nshell_0 = [int(x) for x in other_args[2:]]
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0,
            nmax=nmax_0, n1=n1_0, n2=n2_0, nshell=nshell_0,
            int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    elif len(other_args) == 7:
        a_range0 = list(range(int(other_args[0]), int(other_args[1])+1))
        nmax_0, n1_0, n2_0, nshell_0, ncomponent_0 = [int(x)
                                                      for x in other_args[2:]]
        ncsd_vce_calculations(
            a_prescriptions=a_prescriptions0, a_range=a_range0,
            nmax=nmax_0, n1=n1_0, n2=n2_0,
            nshell=nshell_0, ncomponent=ncomponent_0,
            int_scalefactor=scalefactor0,
            force_ncsd=f_ncsd, force_trdens=f_trdens, force_vce=f_vce,
            force_all=f_all, verbose=verbose0, progress=progress0,
            cluster_submit=cluster_submit0, walltime=walltime0,
        )
    else:
        raise InvalidNumberOfArgumentsException(
            '%d' % (len(argv) - 1,) +
            ' is not a valid number of arguments for ncsm_vce_calc.py.' +
            'Please enter 3-10 arguments.'
        )
