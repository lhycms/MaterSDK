import os
from setuptools import setup
from torch.utils.cpp_extension import BuildExtension, CppExtension


MATERSDK_ROOT = "/data/home/liuhanyu/hyliu/code/matersdk"


se4pw_module = CppExtension(
        name="se4pw_op",
        sources=[
            os.path.join(MATERSDK_ROOT, "source/matersdk/feature/deepmd/src/se4op.cc"),
            os.path.join(MATERSDK_ROOT, "source/matersdk/feature/deepmd/bind/se4pw_op_bind.cc")
        ],
        is_python_module=True,
        extra_compile_args=["-std=c++14"]
)


setup(
    name="se4pw_op",
    ext_modules=[se4pw_module],
    cmdclass={"build_ext": BuildExtension}
)