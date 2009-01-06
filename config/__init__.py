import os

files = os.listdir(os.path.dirname(__file__))
modules = map(lambda m: os.path.splitext(m)[0], 
    filter(lambda f: f.endswith('.py') and not f=='__init__.py', files))
modules.append('ptprofiles')
modules.append('mcasym')
modules.append('away2')
__import__(__name__, {}, {}, modules)

