from propagate import propagate_to

from spice_loader import *

def make_ephemeris(spk_filename, segment_id, *args,
                   type           = 13,
                   append         = False,
                   body_id        = -5440,
                   center_id      = 399,
                   ref_frame_name = 'J2000',
                   degree         = 15,
                   **kwargs):
    if type == 13:
        ts, xs, xf, Phi = propagate_to(*args, **kwargs)
    elif type == 12:
        raise NotImplemented("need function for propagating at a fixed step size or else interpolating a non-fixed step size")

    # Open file
    if append:
        spk = spice.spkopa(spk_filename)
    else: # new file
        if os.path.isfile(spk_filename):
            os.remove(spk_filename)
        spk = spice.spkopn(spk_filename, "SPK_file", 0)


    if len(ts) < degree:
        raise ValueError("polynomial degree is too high for this number of data points")
        
    # Write file
    if type == 13:
        spice.spkw13(spk, body_id, center_id, ref_frame_name, ts[0], ts[-1], segment_id,
                     degree, len(ts), xs.T.copy() / 1000.0, ts.flatten())

    # Close file
    spice.spkcls(spk)

    return ts, xs, xf, Phi
