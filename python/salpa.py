import os
from pathlib import Path
from tempfile import NamedTemporaryFile

_salpa_cmd = Path(__file__).parent.parent.joinpath('build/salpa')


def params(fs_Hz=30000, probechans=384, extrachans=0,
           skip=0, lim=None, thr_rms=3.0, halfwidth=3.0,
           forcepeg=0.3, postblank=0.8,
           postblankinterruptible=False):
    """
    Pre-specify parameters for a SALPA run.

    :param fs_Hz: Sampling rate of recording.
    :type fs_Hz: integer
    :param probechans: Number of electrode channels in recording.
    :type probechans: integer
    :param extrachans: Number of auxiliary channels in recording.
    :type extrachans: integer
    :param skip: Number of scans to skip from beginning of file.
    :type skip: integer
    :param lim: Number of scans to process. (Zero or None means: unlimited.)
    :type lim: integer or None
    :param thr_rms: Threshold for resuming normal operation after an artifact, specified in units of RMS noise for each channel.
    :type thr_rms: float
    :param halfwidth: Half-width of the artifact fitter, in milliseconds.
    :type halfwidth: float
    :param forcepeg: Enforced minimum duration of artifact, in milliseconds.
    :type forcepeg: float
    :param postblank: Number of milliseconds of signal to be blanked after nominal end of artifact.
    :type postblank: float
    :param postblankinterruptible: If given, “blanking” specified by *postblank* is terminated early if the signal on a given channel crosses zero.
    :type postblankinterruptible: bool
    :return: The parameters, combined into a dict that can be passed to ``run()``.
    :rtype: dict

"""
    return {
        'fs_Hz': fs_Hz,
        'probechans': probechans,
        'totalchans': probechans+extrachans,
        'skip': skip,
        'lim': lim,
        'thr_rms': thr_rms,
        'halfwidth': halfwidth,
        'forcepeg': forcepeg,
        'postblank': postblank,
        'postblankinterruptible': postblankinterruptible
    }


def _createpegfile(ss_stim, fs_Hz, forcepeg_ms):
    forcepeg = int(fs_Hz*forcepeg_ms/1000)
    tempfd = NamedTemporaryFile()
    for s in ss_stim:
        tempfd.write(f'{int(s)} {forcepeg}\n')
    return tempfd


def run(ifn, ofn, pars, stim_ss=None):
    """
    Run SALPA on a file

    :param ifn: The file to use for input
    :type ifn: str
    :param ofn: The file to use for output. (Warning: any preexisting file is overwritten.)
    :type ofn: str
    :param pars: Parameters for the SALPA process, constructed with ``params()``.
    :type pars: dict
    :param stim_ss: Optional: list of sample times at which an artifact is to be assumed. (E.g., times of electrical stimuli.) Times are measured in samples, not seconds.
    :type stim_ss: list or np.ndarray
    :return: Nothing
    :exception: If the SALPA program fails to complete, an exception is raised.

"""
    if os.path.exists(ofn):
        os.unlink(ofn)

    if stim_ss is not None:
        stimfd = _createpegfile(stim_ss, pars["fs_Hz"], pars["prepeg"], pars["forcepeg"])

    args = [
        salpa_cmd,
        f'-F {int(pars["fs_Hz"])//1000}',
        f'-c {pars["probechans"]}',
        f'-C {pars["totalchans"]}',
        f'-x {pars["thr_rms"]}',
        f'-l {pars["halfwidth"]}',
        f'-T 12',
        f'-M {pars["skip"]}',
        f'-i "{ifn}"',
        f'-o "{ofn}"',
        f'-f {1}',
        f'-b {pars["postblank"]}'
        ]
    if pars["lim"] is not None:
        args.append(f'-N {pars["lim"]}')
    if stim_ss is not None:
        stimfd = _createpegfile(stim_ss, pars["prepeg"], pars["forcepeg"])
        args.append(f'-P "{stimfd.name}"')
    if not pars["postblankinterruptible"]:
        args.append("-Z")
    if os.system(' '.join(args)):
        print("Failure running SALPA with arguments: ", args)
        raise Exception('SALPA failed')
