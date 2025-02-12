"""
Extracted data processing of DDanalysis
"""
import numpy as np


def apply_data(data_input, correct_lo, crop_w3_min=0, crop_w3_max=0):
    data = {it: data_input[it] for it in ["a1", "a2", "a3"]}  # keep axis the same memory object

    # apply correctLO
    if data_input["LO"] is not None and correct_lo:
        LOi = data_input["LO"].sum(1)
        # note that DDphasing already divides data by LO**0.5 in one direction (likely w3)
        # so v1 makes more sense than v2, although LO ~ E**2 and 2DES ~ E**3
        # also relative peak amplitudes are not corrected for laser spectral profile in w1 direction
        # and detected LO goes through the sample

        # v1
        cor = LOi.reshape((-1, 1, 1)) / LOi.max()
        for it in ["total", "rephasing", "nonrephasing"]:
            if it in data_input:
                data[it] = data_input[it] / cor  # make copy of data, do not modify original data
        data["LO"] = data_input["LO"] / LOi.reshape((-1, 1)) * LOi.max()

        # ~ data = {"total":T/cor, "rephasing":R/cor, "nonrephasing":N/cor, "a1":a1, "a2":a2, "a3":a3, "LO":LO/LOi.reshape((-1,1))*LOi.max()}
        # v2
        # ~ cor = ((LOi/LOi[0])**(3/2)).reshape((-1,1,1))
        # ~ data = {"total":T/cor, "rephasing":R/cor, "nonrephasing":N/cor, "a1":a1, "a2":a2, "a3":a3, "LO":LO/LOi.reshape((-1,1))}
        # v3

        # ~ print("_loadData", T.shape, LO.shape)
        # ~ #T.shape = t2, w1, w3
        # ~ #  w3 is the same as for LO
        # ~ #  but w1 is not, we need to interpolate LO on w1 axis to calculate v2
        # ~ LO1 = np.array([np.interp(a1, a3, it) for it in LO])
        # ~ cor = np.array([np.outer(it1, it2) for it1,it2 in zip(LO1, LO)])
        # ~ data = {"total":T/cor, "rephasing":R/cor, "nonrephasing":N/cor, "a1":a1, "a2":a2, "a3":a3, "LO":LO/LO}

        # v4 from 2023 VCP2 session sand_compensate_lo
        lo_w = data["a3"]
        lo_a = data_input["LO"]

        lo_interpolated_1 = np.asarray([np.interp(data["a1"], lo_w, it) for it in lo_a])
        lo_interpolated_3 = lo_a # this is already interpolated

        # the compensation factor
        COMPENSATION_SCALE = 1
        factor = [np.outer(COMPENSATION_SCALE * lo_1 ** -0.5, (lo_3 / lo_3.max()) ** -0.333)
                     for (lo_3, lo_1) in zip(lo_interpolated_3, lo_interpolated_1)]

        for it in ["total", "rephasing", "nonrephasing"]:
            if it in data_input:
                data[it] = data_input[it] * factor  # make copy of data, do not modify original data
        data["LO"] = data_input["LO"]

    else:
        for it in ["total", "rephasing", "nonrephasing", "LO"]:
            if it in data_input:
                data[it] = data_input[it]  # in this case we can keep the original data, but it must not be modified after
        # ~ data = {"total":T, "rephasing":R, "nonrephasing":N, "a1":a1, "a2":a2, "a3":a3, "LO":LO}

    # apply post crop
    w3n = crop_w3_min
    w3x = crop_w3_max
    if w3n or w3x:  # i.e. non zero
        a3 = data["a3"]
        mn = max(w3n, a3[0]) if (0 < w3n < a3[-1]) else a3[0]
        mx = min(w3x, a3[-1]) if (w3x > 0 and w3x > a3[0]) else a3[-1]
        I = np.logical_and(mn <= a3, a3 <= mx)
        for it in ["total", "rephasing", "nonrephasing", "LO", "a3"]:
            if it in data:
                data[it] = data[it][..., I]
    return data


def amplitude_to_noise(x):
    print("amplitude_to_noise", x.shape)
    N = 2
    sr = np.empty((x.shape[1] - N, x.shape[2] - N), float)
    si = np.empty((x.shape[1] - N, x.shape[2] - N), float)
    S = []
    res = np.abs(x)
    for k in range(x.shape[0]):
        data = x[k]
        for i in range(sr.shape[0]):
            for j in range(sr.shape[1]):
                sr[i, j] = data.real[i:i + N, j:j + N].std()
                si[i, j] = data.imag[i:i + N, j:j + N].std()

        def filter(y, M=1):
            for i in range(M):
                c, edges = np.histogram(y, bins=20, range=(np.nanmin(y), np.nanmax(y)))
                ci = c.argmax()
                my = 0.5 * (edges[ci] + edges[ci + 1])
                sy = np.nanmean(y)
                y[y > my + 1.5 * sy] = np.NaN

        filter(sr, 5)
        filter(si, 5)

        s = (np.nanmean(sr)**2+np.nanmean(si)**2)**0.5
        print("amplitude_to_noise", k, s)
        res[k] /= s
    return res

# WIP - TODO: extract rest of the data processing to here
