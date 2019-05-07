# (c) 2019 Trina Rutz 
# get_sample() inspired by Bart Massey's 'findpeak.py'
# Found at https://github.com/BartMassey/pdx-cs-sound/blob/master/findpeak.py
# goertzel() code helped by Pat Rademacher

import sys
import wave
import math
import numpy as np
import binascii

def write_message(sentance):
    f = open('MESSAGE.txt', 'w')
    f.write(sentance)
    f.close()

# Take every 10 bits: 1st bit is a 0 (start bit, space), end bit is a 1 (stop bit, mark) 
# 8 bits inbetween are an ASCII character
def make_sentance_from_binary(binary):
    sentance_array = []
    for i in range(0, len(binary), 10):
        bin_ascii = []
        if (binary[i] != 0) or (binary[i+9] != 1): 
            raise ValueError('start and stop bit not correct!')
        for j in range(0, 8):
            bin_ascii.append(binary[i + 8 - j])
        letter = chr(int(''.join(str(x) for x in bin_ascii), 2))
        sentance_array.append(letter)
        sentance = ''.join(sentance_array)
    return sentance

def goertzel (sample_rate, samples, targ_freq):
    k = (len(samples) * targ_freq) / sample_rate
    w = ((2 * math.pi) / len(samples)) * k
    cosine_vals = np.cos(w)
    coefficient_vals = 2 * cosine_vals
    s0 = 0
    s1 = 0
    s2 = 0
    for i in range(len(samples)):
        s0 = coefficient_vals * s1 - s2 + samples[i]
        s2 = s1
        s1 = s0
    return np.square(s1) + np.square(s2) - s1 * s2 * coefficient_vals

# every 160 elements of sample array corresponds to a bit 
# divide into chunks of 160 elements and call goertzel on each chunk
def get_binary_from_peaks(sample):
    peak_freqs = []
    for i in range(0, int(len(sample)/160)): 
        chunk = sample[i*160 : i*160 + 160]
        # When each chunk is put thru goertzel, 
        # it'll have a peak at 2225 or 2025Hz
        # Highest peak will show if it's a space or mark bit
        space = goertzel(sample_rate, np.asarray(chunk), 2025)
        mark = goertzel(sample_rate, np.asarray(chunk), 2225)
        peak_freqs.append((lambda s, m: 's' if s > m else 'm')(space, mark))
    # 2225Hz mark corresponds to a 1
    # 2025Hz space corresponds to a 0
    binary = list(map(lambda x: 0 if x == 's' else 1, peak_freqs))
    return binary 

# normalizing sample values to [-1, 1]
# (2* (x - min_x) / (max_x - min_x)) - 1
def normalize(sample):
    maximum = 0.0
    minimum = math.inf
    for i in range(0, len(sample)):
        if sample[i] > maximum:
            maximum = sample[i]
        elif sample[i] < minimum:
            minimum = sample[i]

    new = []
    bottom = maximum - minimum
    for i in range(0, len(sample)):
        top = sample[i] - minimum
        new.append((2 * (top/bottom)) - 1)
    return new

# open wavefile, and get sample
# this function is inspired by Bart Massey's code, more info at top
def get_sample():
    wavefile = wave.open(sys.argv[1], 'rb')
    width = wavefile.getsampwidth()
    print(width)
    wave_bytes = wavefile.readframes(wavefile.getnframes()) 
    # add each frame from file to sample array
    sample = []
    for f in range(0, len(wave_bytes), width):
        current_frame = wave_bytes[f : f+width]
        sample_bytes = current_frame[0 : width]
        sample.append(float(int.from_bytes(sample_bytes, byteorder='little', signed=(width>1))))
    wavefile.close()
    return sample

if __name__ == "__main__":

    sample_rate = 48000
    sample = get_sample()
    normalized = normalize(sample)
    binary = get_binary_from_peaks(normalized)
    sentance = make_sentance_from_binary(binary)
    write_message(sentance)
