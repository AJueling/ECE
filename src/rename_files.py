# this script renames fuiles by replacing str1 with str2
# e.g., called with `python3 rename_files.py /ec/res4/scratch/nkaj/ecearth3/sf02/output/ifs/001 "Sf02" "sf02"`

import os
import sys

if __name__=='__main__':
    dir = sys.argv[1]
    str1 = sys.argv[2]
    str2 = sys.argv[3]
    assert type(str1) is str
    assert type(str2) is str
    for filename in os.listdir(dir):
        # if filename.startswith(str1):
        #     os.rename(f'{dir}/{filename}',f'{dir}/{str2+filename[len(str1):]}')

        if str1 in filename:
            print(filename)
            new_filename = filename.replace(str1,str2)
            print(new_filename)
            os.rename(f'{dir}/{filename}', f'{dir}/{new_filename}')