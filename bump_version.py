#!/usr/bin/env python3
"""
bump_version.py
===============
Bump version number, updating the version.py file in the CWD.
The version.py file should define a __version__ string containing the version
number following the major.minor[.micro] scheme.

Optionally, the version.py can contain a list of files in which
the release "name" should be updated to the new one.
The release name is defined as the release version number preceded
by the letter "v"

e.g. given this version.py file:

__version__ = '1.0.4'
update_files = ['README.md', 'src/main.c']

the actual release name is taken to be 'v1.0.4'.

"""

__author__ = "Simone Marsili"
__version__ = "1.0"

VALID_FIELD_NAMES = ['major', 'minor', 'micro']
NAME_PREFIX = 'v'

def get_command(description=None):
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description,
        epilog=" ")

    parser.add_argument(
        "-f", "--field",
        type=str,
        help="field to bump. valid options are: major minor micro)",
        default="micro")

    return parser.parse_args().field

def bump(string, target_field):
    """Parse version fields in string."""
    import sys

    try:
        try:
            fds = [int(x) for x in string.split('.')]
        except TypeError:
            raise TypeError("Check __version__ (major.minor[.micro] scheme)")
        n_fields = len(fds)
        if n_fields not in [2,3]:
            raise ValueError("Check __version__ (major.minor[.micro] scheme)")
        if n_fields == 2:
            n_fields += 1
            fds.append(0)
        try:
            target_field = VALID_FIELD_NAMES.index(target_field)
        except ValueError:
            raise ValueError("Valid field names are: %s" %
                             ', '.join(VALID_FIELD_NAMES))
        fds = [fds[k]+1 if k == target_field else
               (fds[k] if k < target_field else 0)
               for k in range(3)]
        # if the last field is zero, ignore it
        if fds[-1] == 0:
            fds.pop()
        string = ".".join([str(fd) for fd in fds])
        print("Version bumped to %s" % string, file=sys.stderr)
        return string
    except:
        raise

if __name__ == "__main__":
    import os
    import sys
    field_to_bump = get_command(description=__doc__)
    
    cwd = os.getcwd()
    vpath = cwd + '/version.py'
    # check if version.py is in CWD
    if not os.path.isfile(vpath):
        print("***Type `bump_version.py -h` for help.", file=sys.stderr)
        raise FileNotFoundError(
            "No %s file in %s/: not a valid package." % ('version.py', cwd))

    try:
        import version
    except SyntaxError:
        raise

    try:
        old_version = version.__version__
    except AttributeError:
        # no __version__ attribute in version.py
        raise
    else:
        new_version = bump(old_version, field_to_bump)
        old_name = NAME_PREFIX + old_version
        new_name = NAME_PREFIX + new_version
        # update version file
        with open(vpath, 'w') as f:
            print("__version__ = \'%s\'" % new_version, file=f)

    try:
        update_files = version.update_files
    except AttributeError:
        # no other files were defined in version.py
        update_files = None
    else:
        str_list = '[' + ', '.join(["\'%s\'" % f for f in update_files]) + ']'
        with open(vpath, 'a') as f:
            print("update_files = %s" % str_list, file=f)
        for tfile in update_files:
            print(tfile)
            with open(tfile, 'r') as f:
                tstr = f.read().strip('\n')
            tstr = tstr.replace(old_name, new_name)
            with open(tfile, 'w') as f:
                print(tstr, file=f)
