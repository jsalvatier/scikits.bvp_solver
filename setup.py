# Author: John Salvatier <jsalvati@u.washington.edu>, 2009.
import setuptools


DISTNAME            = 'scikits.bvp_solver'
DESCRIPTION         = "Python package for solving two-point boundary value problems"
LONG_DESCRIPTION    ="""
        bvp_solver is a Python package for solving two-point boundary value problems that wraps
        a slightly modified BVP_SOLVER (see http://cs.stmarys.ca/~muir/BVP_SOLVER_Webpage.shtml).
	Installation indstructions, a tutorial, examples and documentation can be found at 
	http://packages.python.org/scikits.bvp_solver/. If you have questions or suggestions 
        send an e-mail to the mailing list or me.

        To join the mailing list send an e-mail to scikits-bvp_solver+subscribe@googlegroups.com
        """
MAINTAINER          = 'John Salvatier'
MAINTAINER_EMAIL    = "jsalvatier@gmail.com"
URL                 = "http://packages.python.org/scikits.bvp_solver/"
LICENSE             = "BSD"
VERSION             = "1.1"

classifiers =  ['Development Status :: 5 - Production/Stable',
                'Programming Language :: Python',
                'License :: OSI Approved :: BSD License',
                'Intended Audience :: Science/Research',
                'Topic :: Scientific/Engineering',
                'Topic :: Scientific/Engineering :: Mathematics',
                'Operating System :: OS Independent']

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(DISTNAME, parent_package, top_path,
                           namespace_packages = ['scikits'],
                           version = VERSION,
                           maintainer  = MAINTAINER,
                           maintainer_email = MAINTAINER_EMAIL,
                           description = DESCRIPTION,
                           license = LICENSE,
                           url = URL,
                           long_description = LONG_DESCRIPTION)

    config.add_data_files('scikits/__init__.py')
    config.add_extension('bvp_solverf',
                         sources=['scikits/bvp_solver/lib/lampak.f',
				  'scikits/bvp_solver/lib/BVP_LA.f',
				  'scikits/bvp_solver/lib/BVP_M.f90',
                                  'scikits/bvp_solver/lib/BVP_INTERFACE.f90',
				  'scikits/bvp_solver/lib/bvp_interface.pyf'])

    config.add_data_files('scikits/bvp_solver/examples/*.*')
    config.add_data_files('scikits/bvp_solver/lib/BVP_SOLVER_License.txt')
    config.add_data_files('scikits/bvp_solver/lib/BVP_LA_LicenseInfo.txt')
    return config



if __name__ == "__main__":

    from numpy.distutils.core import setup
    setup(configuration=configuration,
        packages = setuptools.find_packages(),
        include_package_data = True,
        platforms = ["any"],
        requires=["numpy"],
        tests_require = ['nose',],
        test_suite='nose.collector',
        zip_safe = True,
        classifiers =classifiers)
