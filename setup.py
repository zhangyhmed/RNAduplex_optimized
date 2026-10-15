import os
import hashlib
import shlex
import subprocess
import sys
from pathlib import Path
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

prefix = Path(os.environ.get('VIENNA_PREFIX', sys.prefix))
header = prefix / 'include/ViennaRNA/duplex.h'
library = prefix / 'lib/libRNA.a'
source_root=Path(__file__).resolve().parent
revision='core-'+hashlib.sha256(b''.join((source_root/name).read_bytes() for name in
    ('native/core.c','native/energy.h','native/python_module.c'))).hexdigest()
native_flags=shlex.split(os.environ.get('RDO_NATIVE_FLAGS','-O3 -flto -march=native -funroll-loops'))+['-fvisibility=hidden']
class BuildNative(build_ext):
    def run(self):
        if not header.is_file() or not library.is_file():
            raise RuntimeError('Set VIENNA_PREFIX to a ViennaRNA 2.7.2 installation with headers and libRNA.a')
        version = subprocess.check_output([str(prefix/'bin/RNAduplex'),'--version'],text=True).strip()
        if version != 'RNAduplex 2.7.2':
            raise RuntimeError(f'Requires ViennaRNA 2.7.2; found {version}')
        super().run()
        destination=Path(self.get_ext_fullpath('rnaduplex_optimized._native')).parent/'_RNAduplex_2_7_2'
        upstream=Path(__file__).resolve().parent/'native/upstream'
        command=[*shlex.split(os.environ.get('CC','gcc')),'-O3','-fopenmp','-DVERSION="2.7.2"',
                 '-I'+str(prefix/'include'),'-I'+str(upstream),
                 str(upstream/'RNAduplex.c'),str(upstream/'RNAduplex_cmdl.c'),str(library),
                 '-L'+str(prefix/'lib'),'-Wl,-rpath,'+str(prefix/'lib'),
                 '-lm','-lstdc++','-lgomp','-lpthread','-o',str(destination)]
        self.spawn(command)


setup(cmdclass={'build_ext':BuildNative},ext_modules=[Extension('rnaduplex_optimized._native',
    sources=['native/core.c','native/python_module.c'],
    depends=['native/energy.h'],
    define_macros=[('RDO_REVISION','"'+revision+'"'),
                   ('RDO_BUILD_FLAGS','"'+' '.join(native_flags)+'"')],
    include_dirs=[str(prefix/'include')],
    extra_objects=[str(library)],
    library_dirs=[str(prefix/'lib')],
    libraries=['m','stdc++','gomp','pthread'],
    runtime_library_dirs=[str(prefix/'lib')],
    extra_compile_args=native_flags,
    extra_link_args=[*native_flags,'-Wl,-Bsymbolic'],
)])
