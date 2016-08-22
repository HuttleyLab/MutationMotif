#!/usr/bin/env python
import doctest, cogent3.util.unit_test as unittest, sys, os
from cogent3.util.misc import app_path

# edited copy of cogent3's alltests

def my_import(name):
    """Imports a module, possibly qualified with periods. Returns the module.
    
    __import__ only imports the top-level module.
    
    Recipe from python documentation at:
    http://www.python.org/doc/2.4/lib/built-in-funcs.html
    """
    mod = __import__(name)
    components = name.split('.')
    for comp in components[1:]:
        mod = getattr(mod, comp)
    return mod

def suite():
    modules_to_test = [
        'test_apps',
        'test_control',
        'test_entropy',
        'test_util',
        'test_heights',
        'test_motif_count',
        'test_complement',
        ]

    alltests = unittest.TestSuite()
    
    for module in modules_to_test:
        if module.endswith('.rst'):
            module = os.path.join(*module.split(".")[:-1]) + ".rst"
            test = doctest.DocFileSuite(module, optionflags=
                doctest.REPORT_ONLY_FIRST_FAILURE |
                doctest.ELLIPSIS)
        else:
            test = unittest.findTestCases(my_import(module))
        alltests.addTest(test)
    return alltests

class BoobyTrappedStream(object):
    def __init__(self, output):
        self.output = output

    def write(self, text):
        self.output.write(text)
        raise RuntimeError("Output not allowed in tests")
        
    def flush(self):
        pass
        
    def isatty(self):
        return False

if __name__ == '__main__':
    if '--debug' in sys.argv:
        s = suite()
        s.debug()
    else:
        orig = sys.stdout
        if '--output-ok' in sys.argv:
            sys.argv.remove('--output-ok')
        else:
            sys.stdout = BoobyTrappedStream(orig)
        try:
            unittest.main(defaultTest='suite', argv=sys.argv)
        finally:
            sys.stdout = orig
