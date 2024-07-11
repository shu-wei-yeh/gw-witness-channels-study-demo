import os
import yaml
import pandas as pd
from gwpy.timeseries import TimeSeries, TimeSeriesDict

# Load configuration from YAML file
with open('config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Extract configuration values
which_channels = config['which_channels']
start_times = config['start_times']
duration = config['duration']
fl = config['fl']
fh = config['fh']
fftlength = config['fftlength']
overlap = config['overlap']
strain_channel = config['strain_channel']
wit_channels = config['wit_channels']
threshold = config['threshold']
folder_path = config['folder_path']
output_fig_path = config['output_fig_path']

# Read channels from file
with open(wit_channels, 'r') as file:
    channels = [line.strip() for line in file]

# Ensure output directory exists
os.makedirs(output_fig_path, exist_ok=True)

def folder_files(j):
    return f'{folder_path}{int((j - j % 100000) / 100000)}/K-K1_R-{j}-32.gwf'

# Read data for each start time
def read_data(start_time):
    end_time = start_time + duration
    first_file = int(start_time - start_time % 32)
    last_file = int(end_time - end_time % 32)
    
    strain_channel_files = [folder_files(j) for j in range(first_file, last_file + 1, 32)]
    
    strain_data = TimeSeriesDict.read(
        strain_channel_files, 
        strain_channel, 
        start=start_time, 
        end=end_time
    )[strain_channel[0]]
    
    witness_data_dict = TimeSeriesDict.read(
        strain_channel_files, 
        channels, 
        start=start_time, 
        end=end_time
    )
    
    return strain_data, witness_data_dict

def max_coh_calculation(start_times, fl, fh, fftlength, overlap, threshold):
    results = []
    for start_time in start_times:
        strain_data, witness_data_dict = read_data(start_time)
        for channel in channels:
            witness_data = witness_data_dict[channel]
            
            if strain_data.sample_rate != witness_data.sample_rate:
                strain_data_resampled = strain_data.resample(witness_data.sample_rate)
            else:
                strain_data_resampled = strain_data

            coherence = TimeSeries.coherence(
                strain_data_resampled,
                witness_data,
                fftlength=fftlength,
                overlap=overlap,
                window='hann'
            )
            max_coh_value = coherence.crop(fl, fh).max().value
            if max_coh_value >= threshold:
                max_coh = round(max_coh_value, 2)
                results.append([start_time, duration, channel, max_coh])
        results.append([])  # Add a newline for each new start time
    return results

# Calculate and save the maximum coherence values
coh_results = max_coh_calculation(start_times, fl, fh, fftlength, overlap, threshold)

# Save results to a text file
output_file_path = os.path.join(output_fig_path, f'max_coherence_values_of_{which_channels}_channels_{fl}Hz-{fh}Hz-{start_times[0]}.txt')
df = pd.DataFrame(coh_results, columns=['Start Time', 'Duration', 'Channel', 'Max Coherence'])
df.to_csv(output_file_path, sep='\t', index=False, header=True, mode='w', line_terminator='\n')

print(f"Maximum coherence values saved to {output_file_path}")