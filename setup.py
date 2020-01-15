#!/usr/bin/env python3

from setuptools import setup

MAJOR = 0
MINOR = 0
TINY  = 1
version='%d.%d.%d' % (MAJOR, MINOR, TINY)


setup(name='trajectory',
      version=version,
      description='Trajectory analysis library for lunar missions',
      author='John O. Woods, Ph.D.',
      author_email='john@openlunar.org',
      url='http://www.openlunar.org',
      include_package_data=True,
      packages=['test'],
      test_suite='test.trajectory_test_suite'
      )
