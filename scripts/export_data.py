#!/usr/bin/env python
# -*- coding: utf-8 -*-

from jetbreakup import *
import os
import shutil

output_filename  = 'pipe-jet-breakup-data'
output_directory = './'+output_filename

with open(output_filename+'.pickle') as f:
   df_jet_breakup, metadata = pickle.load(f)

# TODO: Add README file with more information about the archive.
# TODO: Change filename of zip file to include today's date. Also add a file to the zip file with SVN repository information, including when last updated and when the zip file was created.
# TODO: Add example files for reading and plotting the data with both Python (CSV and Pandas) and also Matlab (CSV).

if os.path.exists(output_directory):
   shutil.rmtree(output_directory)

os.makedirs(output_directory)
os.makedirs(output_directory+'/photos/')

print 'Writing data frame to CSV file...'
df_jet_breakup.to_csv(output_directory+'/'+output_filename+'.csv', sep='\t', encoding='utf-8')
print 'Copying data frame...'
shutil.copyfile(output_filename+'.pickle', output_directory+'/'+output_filename+'.pickle')
print 'Copying bibliography...'
shutil.copyfile('../pipe-jet-breakup-data.bib', output_directory+'/pipe-jet-breakup-data.bib')

df_photos = df_jet_breakup[df_jet_breakup['photo filename'].notnull()]
photo_file_location_array = []
photo_filename_array      = []
for photo_filename in df_photos['photo filename']:
   item_df = df_photos[df_photos['photo filename'] == photo_filename]
   for key in item_df['key']:
      #print key, photo_filename
      
      photo_file_location = '../data/'+key+'/'+photo_filename
      
      if not(photo_filename in photo_filename_array):
         photo_filename_array.append(photo_filename)
         photo_file_location_array.append(photo_file_location)

i = 0
for photo_filename in photo_filename_array:
   photo_file_location = photo_file_location_array[i]
   
   print 'Copying', photo_filename, '...'
   shutil.copyfile(photo_file_location, output_directory+'/photos/'+photo_filename)

print 'Zipping...'
shutil.make_archive(output_filename, 'zip', output_directory)

print 'Cleaning up...'
shutil.rmtree(output_directory)
