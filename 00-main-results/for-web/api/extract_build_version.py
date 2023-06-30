import pickle
import os

from snapanalysis.config import OUTPUT_DIRECTORY as SNAPANALYSIS_OUTPUT_DIRECTORY

OUTPUT_FILE = os.path.join(SNAPANALYSIS_OUTPUT_DIRECTORY, 'web/api/precompiled/snapanalysis_build.pickle')
VERSION_FILE = '/build/info.txt'

def main():

    with open(VERSION_FILE, 'r') as f:
        version_info = f.read()

    date, commit, __ = version_info.split('\n')

    date = date.partition('=')[-1]
    commit = commit.partition('=')[-1]

    version = {
        'date': date,
        'commit': commit,
    }

    with open(OUTPUT_FILE, 'wb') as f:
        pickle.dump(version, f, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()