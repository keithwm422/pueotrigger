import numpy as np, os, noise
import tools.CoRaLs_geometry as corals_geometry
import coherent_sum

def stream_power_1s(window=160, step=40, target_seconds=10.0,
                    block_fbins=2**22, profile='v2', vrms=1.0,
                    out_file='noise/power_160_40_1s.npy'):
    fs = 1.0 / (corals_geometry.ritc_sample_step * 1e-9)
    total_samples = int(target_seconds * fs)
    os.makedirs('noise', exist_ok=True)


    tn = noise.ThermalNoise(0.26, 0.95, filter_order=(10,10), v_rms=vrms,
                                fbins=block_fbins,
                                time_domain_sampling_rate=corals_geometry.ritc_sample_step)
    
    tail = None
    powers = []
    written = 0
    blk = 0
    while written < total_samples:
        blk += 1
        _, _, wf = tn.makeNoiseWaveform(ntraces=1)
        wf = wf[0].real  # 1D
        take = min(block_fbins, total_samples - written)
        wf = wf[:take]
        written += take

        if tail is not None:
            wf = np.concatenate([tail, wf])
        if len(wf) >= window:
            p, _ = coherent_sum.powerSum(wf, window=window, step=step)
            # Keep last (window + step*2) samples as tail for boundary continuity
            keep = window + step*2
            tail = wf[-keep:] if keep < len(wf) else wf
            powers.append(p.astype('float32'))
        else:
            tail = wf
        print(f"  Block {blk}: samples_accum={written}/{total_samples}  power_frames_total={sum(len(a) for a in powers)}")

    power_all = np.concatenate(powers)
    np.save(out_file, power_all)
    print(f"[done] Saved {out_file} frames={len(power_all)} (~{len(power_all)*step/fs:.3f} s equivalent)")

if __name__ == '__main__':
    stream_power_1s()