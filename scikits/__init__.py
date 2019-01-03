# Activate namespace packages.
__import__('pkg_resources').declare_namespace(__name__)

# Add path to dll library to PATH on Windows
import os, sys
if sys.platform == 'win32':
    dllpath = os.path.join(os.path.dirname(__file__), '.libs')
    os.environ['PATH'] += dllpath + ';'
