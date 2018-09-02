from setuptools import setup
import os
import sys


if sys.version_info[0]!=3:
    print("pykbi requires python 3.")
    sys.exit(1)


def get_requirements():
    """
    Read requirements.txt
    """

    thispath = os.path.abspath(os.path.dirname(__file__))

    requirements = []
    fname = os.path.join(thispath,"requirements.txt")
    with open(fname, encoding='utf-8') as f:
        for line in f:
            requirements.append(line.strip())
    return requirements


setup(
        name="pykbi",
        version="1.0.0",
        description="Integrate Radial distribution functions, and calculate Kirkwood-Buff integrals from open and closed systems",
        keywords="kirkwood-buff integral",
        install_requires=get_requirements(),
        author="Sondre K. Schnell",
        author_email="sondresc@gmail.com",
        packages=["pykbi"],
        )
