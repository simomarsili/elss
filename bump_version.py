#!/usr/bin/env python3
DESCRIPTION =\
"""
bump_version.py
===============
Bump version number. A __version__ variable should be defined in a
version.py file within the package dir.
The version numbering should follow the major.minor[.micro] scheme.
"""

__author__ = "Simone Marsili"
__version__ = "0.1.0"

def get_command(description):
    """syntax: bump_version.py -p <package_dir> [-f major|minor|micro]."""
    import os
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=description,
        epilog=" ")
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument(
        "-p", "--package_dir", type=str,
        help="package directory (type . for present work directory)",
        required=True)
    parser.add_argument(
        "-f", "--field", type=str,
        help="version field to bump. valid options are: major minor micro)",
        default="micro")

    args = parser.parse_args()
    pdir = args.package_dir
    field = args.field

    # chk args
    if not os.path.isdir(pdir):
        raise FileNotFoundError("No such package/directory: %s" % pdir)
    else:
        filename = pdir + "/" + "version.py"
        if not os.path.isfile(filename):
            raise FileNotFoundError(
                "No %s file in %s/: not a valid package." % (filename, pdir))
    return (pdir, field)

def bump(string, target_field):
    """Parse version fields in string."""
    import re
    import sys
    field_values = ["major", "minor", "micro"]
    try:
        VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
        mo = re.search(VSRE, string, re.M)
        try:
            fds = [int(x) for x in mo.group(1).split(".")]
        except TypeError:
            raise TypeError("Version numbering should follow the major.minor[.micro] convention")
        n_fields = len(fds)
        if n_fields not in [2,3]:
            raise ValueError("Version numbering should follow the major.minor[.micro] convention")
        if n_fields == 2:
            n_fields += 1
            fds.append(0)
        try:
            target_field = field_values.index(target_field)
        except ValueError:
            raise ValueError("Possible fields are major/minor/micro")
        fds = [fds[k]+1 if k == target_field else
               (fds[k] if k < target_field else 0)
               for k in range(3)]
        # if the last field is zero, it wont be printed
        if fds[-1] == 0:
            fds.pop()
        string = ".".join([str(fd) for fd in fds])
        print("Version bumped to %s" % string, file=sys.stderr)
        return '__version__ = \'' + string + '\''
    except:
        raise
        return string

if __name__ == "__main__":
    package_dir, field_to_increase = get_command(DESCRIPTION)
    vpath = package_dir + "/version.py"

    with open(vpath, "r") as f:
        lines = f.readlines()
    with open(vpath, "w") as f:
        for line in lines:
            new_version = bump(line, field_to_increase)
            print(new_version, file=f)
