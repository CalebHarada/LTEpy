from setuptools import setup, find_packages
import re

# auto-updating version code stolen from orbitize! stolen from RadVel
def get_property(prop, project):
    result = re.search(
        r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
        open(project + "/__init__.py").read(),
    )
    return result.group(1)

# stolen from orbitize!
def get_requires():
    reqs = []
    for line in open("requirements.txt", "r").readlines():
        reqs.append(line)
    return reqs

setup(
    name="LTEpy",
    version=get_property("__version__", "LTEpy"),
    description="LTEpy is a package for calculating and vizualing properties of a gas in local thermodynamic equilibrium.",
    packages=find_packages(),
    author="Gardiner, E. C., Harada, C.",
    author_email="ecg@berkeley.edu, charada@berkeley.edu",
    license="MIT",
    url="https://github.com/CalebHarada/LTEpy"

)