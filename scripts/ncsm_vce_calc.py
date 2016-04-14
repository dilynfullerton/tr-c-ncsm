#!/usr/bin/python
"""ncsm_vce_calc.py

To run as a script:

    $ ncsm_vce_calc.py [-f[ntv]] [-v] [-s scalefactor] [-t walltime]
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
        force recalculation of NCSD; this also forces TRDENS and VCE
    t
        force recalculation of TRDENS; this also forces VCE
    v
        force recalculation of VCE
-m or -M or -e
    The first two arguments are used to determine a range of A prescriptions
    -m
        Prescriptions are all increasing length-3 combinations
          of integers in the range [Ap_min, Ap_max]
    -M
        Prescriptions are all increasing length-3 combinations with repetition
          of integers in the range [Ap_min, Ap_max]
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
from Queue import Queue
from math import floor
from os import path, remove, link
from subprocess import Popen, PIPE
from sys import argv, stdout
from threading import Thread, currentThread

from FdoVCE import run as vce_calculation
from InvalidNumberOfArgumentsException import InvalidNumberOfArgumentsException

from _make_dirs import DPATH_TEMPLATES, DPATH_RESULTS
from _make_dirs import make_vce_directories, make_trdens_file, rename_egv_file
from _make_dirs import prepare_directories, remove_ncsd_tmp_files
from _make_dirs import EgvFileNotFoundException


# CONSTANTS
N_SHELL = 1
N_COMPONENT = 2
NMAX = 6
N1 = 15
N2 = 15

# files
FNAME_FMT_VCE_INT = 'A%d.int'  # A value
FNAME_TRDENS_OUT = 'trdens.out'
FNAME_HEFF = 'Heff_OLS.dat'

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

# other
WIDTH_TERM = 79
WIDTH_PROGRESS_BAR = 48
STR_PROGRESS_BAR = '    Progress: %3d/%-3d '
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
    # Check that files exist
    for f in [he4_fpath, he5_fpath, he6_fpath]:
        if not path.exists(f):
            raise NcsdOutfileNotFoundException(
                'NCSD outfile not found: %s' % f)
    # Do vce
    a_0 = a_range[0]
    fpath_fmt = path.join(dirpath_vce, FNAME_FMT_VCE_INT)
    fpath = fpath_fmt % a_0
    if force or not path.exists(fpath):
        vce_calculation(a_prescription, fpath, he4_fpath, he5_fpath, he6_fpath)
    # Same interaction is used for all masses
    if len(a_range) > 1:
        for a in a_range[1:]:
            next_fpath = fpath_fmt % a
            if path.exists(next_fpath) and force:
                remove(next_fpath)
                link(fpath, next_fpath)
            elif not path.exists(next_fpath):
                link(fpath, next_fpath)


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


def _ncsd_multiple_calculations_t(
        a_aeff_set, a_aeff_to_dpath_map, a_aeff_to_egvfile_map,
        force, progress=True, str_prog_ncsd=STR_PROG_NCSD,
        max_open_threads=MAX_OPEN_THREADS
):
    def _r(a_, aeff_, q):
        _run_ncsd(
            dpath=a_aeff_to_dpath_map[(a_, aeff_)],
            fpath_egv=a_aeff_to_egvfile_map[(a_, aeff_)],
            force=force, verbose=False
        )
        q.put(currentThread())

    todo_list = list(a_aeff_set)
    active_list = list()
    done_queue = Queue()

    jobs_completed = 0
    jobs_total = len(todo_list)

    completed_dpath_list = list()
    thread_dpath_map = dict()

    if progress and jobs_total > 0:
        print str_prog_ncsd
        _print_progress(jobs_completed, jobs_total)
    while len(todo_list) > 0 or len(active_list) > 0:
        # if room, start new threads
        while len(todo_list) > 0 and len(active_list) < max_open_threads:
            a, aeff = todo_list.pop()
            t = Thread(target=_r, args=(a, aeff, done_queue))
            active_list.append(t)
            thread_dpath_map[t] = a_aeff_to_dpath_map(a, aeff)
            t.start()
        # remove any threads that have finished
        if not done_queue.empty():
            while not done_queue.empty():
                t = done_queue.get()
                t.join()
                active_list.remove(t)
                completed_dpath_list.append(thread_dpath_map[t])
                jobs_completed += 1
            if progress:
                _print_progress(jobs_completed, jobs_total)
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)
    return completed_dpath_list


def _ncsd_multiple_calculations_s(
        a_aeff_set,
        a_aeff_to_dpath_map, a_aeff_to_egvfile_map, a_aeff_to_jobfile_map,
        force, progress,
        fname_stdout=FNAME_QSUB_STDOUT, fname_stderr=FNAME_QSUB_STDERR,
):
    submitted_jobs = 0
    if progress and len(a_aeff_set) > 0:
        print '  Submitting jobs...'
    completed_dpath_list = list()
    for a, aeff in a_aeff_set:
        fpath_egv = a_aeff_to_egvfile_map[(a, aeff)]
        if path.exists(fpath_egv) and not force:
            completed_dpath_list.append(a_aeff_to_dpath_map[(a, aeff)])
        else:
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
    if progress:
        if submitted_jobs == 1:
            print '  1 job submitted to cluster'
        else:
            print '  %d jobs submitted to cluster' % submitted_jobs
    return completed_dpath_list


def _ncsd_multiple_calculations(
        a_aeff_set, a_aeff_to_dpath_map, a_aeff_to_egvfile_map,
        force, verbose, progress, str_prog_ncsd=STR_PROG_NCSD
):
    # do ncsd
    jobs_total = len(a_aeff_set)
    jobs_completed = 0
    progress = progress and not verbose
    if progress and jobs_total > 0:
        print str_prog_ncsd
    completed_dpath_list = list()
    for a, aeff in sorted(a_aeff_set):
        dpath = a_aeff_to_dpath_map[(a, aeff)]
        if progress:
            _print_progress(jobs_completed, jobs_total)
        _run_ncsd(
            dpath=dpath, fpath_egv=a_aeff_to_egvfile_map[(a, aeff)],
            force=force, verbose=verbose
        )
        completed_dpath_list.append(dpath)
        jobs_completed += 1
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)
    return completed_dpath_list


def ncsd_multiple_calculations(
        a_presc_list, a_values, z, nmax, a_0,
        nshell=N_SHELL,  n1=N1, n2=N1, scalefactor=None,
        force=False, verbose=False, progress=True, threading=True,
        cluster_submit=False, walltime=None, remove_tmp_files=True,
        str_prog_ncsd=STR_PROG_NCSD,
        dpath_templates=DPATH_TEMPLATES, dpath_results=DPATH_RESULTS,
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
    :param remove_tmp_files: if true, removes all of the remnant *.tmp files
    following the NCSD calculation
    :param str_prog_ncsd: string to show before progress bar
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    """
    # make (A, Aeff, Nhw) set
    a_aeff_nhw_set = set()
    for ap in a_presc_list:
        a_aeff_nhw_set |= set(zip(
            a_values, ap, [nmax + a - a_0 for a in a_values]
        ))
    # separate set into lists
    a_list, aeff_list, nhw_list = list(), list(), list(),
    for a, aeff, nhw in a_aeff_nhw_set:
        a_list.append(a)
        aeff_list.append(aeff)
        nhw_list.append(nhw)
    # prepare directories and get maps
    a_aeff_maps = prepare_directories(
        a_list=a_list, aeff_list=aeff_list, nhw_list=nhw_list,
        z=z, n1=n1, n2=n2, nshell=nshell, scalefactor=scalefactor,
        cluster_submit=cluster_submit, walltime=walltime, progress=progress,
        dpath_results=dpath_results, dpath_templates=dpath_templates,
        force=force,
    )
    a_aeff_to_dir, a_aeff_to_egv, a_aeff_to_job, a_aeff_to_out = a_aeff_maps
    # make (A, Aeff) set and do NCSD
    a_aeff_set = set([(a, aeff) for a, aeff, nhw in a_aeff_nhw_set])
    if cluster_submit:
        completed_dpath_list = _ncsd_multiple_calculations_s(
            a_aeff_set=a_aeff_set,
            a_aeff_to_dpath_map=a_aeff_to_dir,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            a_aeff_to_jobfile_map=a_aeff_to_job,
            progress=progress, force=force,
        )
    elif threading and len(a_aeff_nhw_set) > 1:
        completed_dpath_list = _ncsd_multiple_calculations_t(
            a_aeff_set=a_aeff_set,
            a_aeff_to_dpath_map=a_aeff_to_dir,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            force=force, progress=progress,
            str_prog_ncsd=str_prog_ncsd
        )
    else:
        completed_dpath_list = _ncsd_multiple_calculations(
            a_aeff_set=a_aeff_set,
            a_aeff_to_dpath_map=a_aeff_to_dir,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            force=force, verbose=verbose, progress=progress,
            str_prog_ncsd=str_prog_ncsd
        )
    if remove_tmp_files:
        remove_ncsd_tmp_files(completed_dpath_list)
    return a_aeff_maps


def ncsd_single_calculation(
        a, aeff, z, scalefactor,
        nhw=NMAX, n1=N1, n2=N1, nshell=N_SHELL,
        force=False, verbose=False, progress=True,
        cluster_submit=False, walltime=None, remove_tmp_files=True,
        dpath_templates=DPATH_TEMPLATES, dpath_results=DPATH_RESULTS,
):
    """Performs a single NCSD calculation for the given (A, Aeff) pair
    :param a: actual mass number A
    :param aeff: effective mass number Aeff
    :param z: proton number Z
    :param scalefactor: factor by which to scale off-diagonal coupling terms of
    the TBME interaction
    :param nhw: model space truncation
    :param n1: max allowed 1-particle state
    :param n2: max allowed 2-particle state
    :param nshell: shell (0: s, 1: p, 2: sd, ...)
    :param force: if true, does the calculation even if output files already
    exist
    :param verbose: if true, prints the regular output of NCSD to stdout, else
    suppresses this and saves it in a file
    :param progress: if true, shows a progress bar (verbose mode will be off)
    :param cluster_submit: if true, submits the job to the OpenMP cluster
    using qsub
    :param walltime: walltime to be allotted to a cluster submission
    :param remove_tmp_files: if true, removes remnant *.tmp files following
    the NCSD calculation (Note this will not be the case for a job submitted
    to the cluster)
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    """
    a_aeff_to_dpath, a_aeff_to_egv, a_aeff_to_job = prepare_directories(
        a_list=[a], aeff_list=[aeff], nhw_list=[nhw],
        z=z, n1=n1, n2=n2, nshell=nshell, scalefactor=scalefactor,
        cluster_submit=cluster_submit, walltime=walltime, progress=progress,
        dpath_templates=dpath_templates, dpath_results=dpath_results,
        force=force,
    )[:3]
    if cluster_submit:
        completed_dpaths_list = _ncsd_multiple_calculations_s(
            a_aeff_set=set([(a, aeff)]),
            a_aeff_to_dpath_map=a_aeff_to_dpath,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            a_aeff_to_jobfile_map=a_aeff_to_job,
            force=force, progress=progress,
        )
    else:
        completed_dpaths_list = _ncsd_multiple_calculations(
            a_aeff_set=set([(a, aeff)]),
            a_aeff_to_dpath_map=a_aeff_to_dpath,
            a_aeff_to_egvfile_map=a_aeff_to_egv,
            force=force, progress=progress,
            verbose=verbose,
        )
    if remove_tmp_files:
        remove_ncsd_tmp_files(completed_dpaths_list)


def ncsd_exact_calculations(
        z, a_range,
        nmax=NMAX, n1=N1, n2=N2,
        nshell=N_SHELL, ncomponent=N_COMPONENT, int_scalefactor=None,
        force=False, verbose=False, progress=True,
        cluster_submit=False, walltime=None, remove_tmp_files=True,
        str_prog_ncsd_ex=STR_PROG_NCSD_EX,
        dpath_templates=DPATH_TEMPLATES, dpath_results=DPATH_RESULTS,
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
    :param remove_tmp_files: if true, removes all of the remnant *.tmp files
    following the NCSD calculation
    :param str_prog_ncsd_ex: string to show before progress bar
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    """
    ncsd_multiple_calculations(
        z=z, a_values=a_range, a_presc_list=[a_range],
        nmax=nmax, n1=n1, n2=n2, nshell=nshell, scalefactor=int_scalefactor,
        a_0=_generating_a_values(n_shell=nshell, n_component=ncomponent)[0],
        cluster_submit=cluster_submit, walltime=walltime,
        force=force, verbose=verbose, progress=progress,
        remove_tmp_files=remove_tmp_files,
        str_prog_ncsd=str_prog_ncsd_ex,
        dpath_templates=dpath_templates, dpath_results=dpath_results,
    )


class NcsdOutfileNotFoundException(Exception):
    pass


def vce_single_calculation(
        a_values, a_prescription, a_range, z, nmax,
        a_aeff_dir_map, a_aeff_outfile_map,
        n1=N1, n2=N1, nshell=-1, ncomponent=-1, int_scalefactor=None,
        force_trdens=False, force_vce=False, verbose=False,
        dpath_results=DPATH_RESULTS, dpath_templates=DPATH_TEMPLATES,
):
    """Valence cluster expansion
    :param a_values: 3-tuple of A values that form the base for constructing
    Heff
    :param a_prescription: 3-tuple of Aeff values used in place of the actual
    A values in constructing Heff
    :param a_range: sequence of A values for which effective interaction files
    are to be generated
    :param z: protonn number
    :param nmax: major oscillator model space truncation
    :param a_aeff_dir_map: map from (A, Aeff) tuple to the directory in which
    this calculation is being done
    :param a_aeff_outfile_map: map from (A, Aeff) tuple to the *.out file
    produced by NCSD for this pair
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
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    """
    # check that files exist
    for f in a_aeff_outfile_map.values():
        if not path.exists(f):
            raise NcsdOutfileNotFoundException(
                'NCSD outfile not found: %s' % f)
    # for the 3rd a value, make trdens file and run TRDENS
    a_aeff6 = (a_values[2], a_prescription[2])
    make_trdens_file(z=z, a=a_values[2], nuc_dir=a_aeff_dir_map[a_aeff6],
                     dpath_results=dpath_results, dpath_temp=dpath_templates)
    a6_dirpath = a_aeff_dir_map[a_aeff6]
    try:
        rename_egv_file(
            a6_dir=a6_dirpath, nhw=nmax+2, a6=a_values[2], force=force_trdens)
    except EgvFileNotFoundException:
        raise
    _run_trdens(a6_dir=a6_dirpath, force=force_trdens, verbose=verbose)

    # do valence cluster expansion
    vce_dirpath = make_vce_directories(
        a_prescription=a_prescription, nmax=nmax, n1=n1, n2=n2,
        nshell=nshell, ncomponent=ncomponent, scalefactor=int_scalefactor,
        dpath_results=dpath_results,
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
        a_aeff_to_dpath_map, a_aeff_to_out_fpath_map,
        dpath_templates, dpath_results,
        force_trdens, force_vce, verbose, progress,
        int_scalefactor=None, max_open_threads=MAX_OPEN_THREADS,
        str_prog_vce=STR_PROG_VCE,
):
    error_messages = Queue()

    def _r(ap0, q):
        try:
            vce_single_calculation(
                z=z, a_values=a_values,
                a_prescription=ap0, a_range=a_range,
                nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
                int_scalefactor=int_scalefactor,
                force_trdens=force_trdens, force_vce=force_vce, verbose=verbose,
                dpath_templates=dpath_templates, dpath_results=dpath_results,
                a_aeff_dir_map=a_aeff_to_dpath_map,
                a_aeff_outfile_map=a_aeff_to_out_fpath_map,
            )
            q.put(currentThread())
        except EgvFileNotFoundException, e:
            error_messages.put(str(e))
            q.put(currentThread())
        except NcsdOutfileNotFoundException, e:
            error_messages.put(str(e))
            q.put(currentThread())

    todo_list = list(a_presc_list)
    active_list = list()
    done_queue = Queue()

    jobs_completed = 0
    jobs_total = len(todo_list)

    if progress and jobs_total > 0:
        print str_prog_vce
        _print_progress(jobs_completed, jobs_total)
    while len(todo_list) > 0 or len(active_list) > 0:
        # if room, start new threads
        while len(todo_list) > 0 and len(active_list) < max_open_threads:
            ap = todo_list.pop()
            t = Thread(target=_r, args=(ap, done_queue))
            active_list.append(t)
            t.start()
        # remove any threads that have finished
        if not done_queue.empty():
            while not done_queue.empty():
                t = done_queue.get()
                t.join()
                active_list.remove(t)
                jobs_completed += 1
            if progress:
                _print_progress(jobs_completed, jobs_total)
    if progress:
        _print_progress(jobs_completed, jobs_total, end=True)
    while not error_messages.empty():
        em = error_messages.get()
        print em
    return jobs_total == jobs_completed


def _vce_multiple_calculations(
        z, a_values, a_presc_list, a_range, nmax, n1, n2, nshell, ncomponent,
        a_aeff_to_out_fpath_map, a_aeff_to_dpath_map,
        dpath_templates, dpath_results,
        force_trdens, force_vce, verbose, progress,
        int_scalefactor=None,
        str_prog_vce=STR_PROG_VCE,
):
    jobs_total = len(a_presc_list)
    jobs_completed = 0
    progress = progress and not verbose
    if progress and jobs_total > 0:
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
                dpath_results=dpath_results, dpath_templates=dpath_templates,
                a_aeff_outfile_map=a_aeff_to_out_fpath_map,
                a_aeff_dir_map=a_aeff_to_dpath_map,
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
    return jobs_completed == jobs_total


def vce_multiple_calculations(
        a_values, a_presc_list, a_range, z, nmax, n1, n2, nshell, ncomponent,
        a_aeff_to_out_fpath_map, a_aeff_to_dpath_map,
        force_trdens, force_vce, verbose, progress, threading,
        int_scalefactor=None,
        dpath_templates=DPATH_TEMPLATES, dpath_results=DPATH_RESULTS,
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
    :param a_aeff_to_dpath_map: map from (A, Aeff) to the directory in which
    this calculation is being done
    :param a_aeff_to_out_fpath_map: map from (A, Aeff) to the *.out file
    produced by NCSD for this pair
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
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    """
    if threading and len(a_presc_list) > 1:
        _vce_multiple_calculations_t(
            a_values=a_values, a_presc_list=a_presc_list, a_range=a_range,
            z=z, nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
            int_scalefactor=int_scalefactor,
            force_trdens=force_trdens, force_vce=force_vce,
            verbose=verbose, progress=progress,
            dpath_results=dpath_results, dpath_templates=dpath_templates,
            a_aeff_to_out_fpath_map=a_aeff_to_out_fpath_map,
            a_aeff_to_dpath_map=a_aeff_to_dpath_map,
        )
    else:
        return _vce_multiple_calculations(
            a_values=a_values, a_presc_list=a_presc_list, a_range=a_range,
            z=z, nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
            int_scalefactor=int_scalefactor,
            force_trdens=force_trdens, force_vce=force_vce,
            verbose=verbose, progress=progress,
            dpath_results=dpath_results, dpath_templates=dpath_templates,
            a_aeff_to_out_fpath_map=a_aeff_to_out_fpath_map,
            a_aeff_to_dpath_map=a_aeff_to_dpath_map,
        )


def ncsd_vce_calculations(
        a_prescriptions, a_range,
        nmax=NMAX, n1=N1, n2=N2, nshell=N_SHELL, ncomponent=N_COMPONENT,
        int_scalefactor=None,
        force_ncsd=False, force_trdens=False, force_vce=False, force_all=False,
        verbose=False, progress=True, threading=True,
        cluster_submit=False, walltime=None, remove_tmp_files=True,
        dpath_templates=DPATH_TEMPLATES, dpath_results=DPATH_RESULTS,
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
    :param remove_tmp_files: if true, removes all of the *.tmp files following
    the ncsd calculations
    :param dpath_templates: path to the templates directory
    :param dpath_results: path to the results directory
    """
    a_values = _generating_a_values(n_shell=nshell, n_component=ncomponent)
    z = int(a_values[0] / ncomponent)
    a_presc_list = list(a_prescriptions)
    a_aeff_maps = ncsd_multiple_calculations(
        z=z, a_values=a_values, a_presc_list=a_presc_list,
        nmax=nmax, a_0=a_values[0], n1=n1, n2=n2, nshell=nshell,
        scalefactor=int_scalefactor,
        force=force_all or force_ncsd,
        verbose=verbose, progress=progress, threading=threading,
        cluster_submit=cluster_submit, walltime=walltime,
        remove_tmp_files=remove_tmp_files,
        dpath_templates=dpath_templates, dpath_results=dpath_results,
    )
    a_aeff_to_dir, a_aeff_to_egv, a_aeff_to_job, a_aeff_to_out = a_aeff_maps
    vce_multiple_calculations(
        z=z, a_values=a_values, a_presc_list=a_presc_list, a_range=a_range,
        nmax=nmax, n1=n1, n2=n2, nshell=nshell, ncomponent=ncomponent,
        int_scalefactor=int_scalefactor,
        force_trdens=force_trdens or force_all,
        force_vce=force_vce or force_all,
        verbose=verbose, progress=progress, threading=threading,
        dpath_results=dpath_results, dpath_templates=dpath_templates,
        a_aeff_to_out_fpath_map=a_aeff_to_out,
        a_aeff_to_dpath_map=a_aeff_to_dir,
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
        force_trdens = True
        force_vce = True
    if 't' in argv0:
        force_trdens = True
        force_vce = True
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
