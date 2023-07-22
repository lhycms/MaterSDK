from matersdk.io.publicLayer.structure import DStructure

structure = DStructure.from_file(
                file_format="pwmat", 
                file_path="/data/home/liuhanyu/hyliu/code/matersdk/demo/structure/atom.config")
print(structure.cart_coords)