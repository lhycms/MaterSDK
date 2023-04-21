#!/Users/mac/opt/anaconda3/envs/research/bin/python3
import click
from matersdk.io.publicLayer.structure import DStructure


@click.command()
@click.option('--input_file',
            #default=1,
            type=str,
            prompt='The path of input file',
            help="The path of input file.",
            )
@click.option('--input_format',
            type=str,
            prompt='The format of input file',
            help='The format of input file.',
            )
@click.option("--output_file",
            type=str,
            prompt="The path of output file",
            help="The path of output file.",
            )
@click.option('--output_format',
            type=str,
            prompt='The format of output file',
            help='The format of output file.',
            )
def convert_format(input_file, input_format, output_file, output_format):
    """
    Convert the format of `structure file` between: \n
        1. pwmat
        2. poscar
        3. cssr
        4. json
        5. xsf
        6. mcsqs
        7. prismatic
        8. yaml
        9. fleur-inpgen
    """
    structure = DStructure.from_file(file_path=input_file, file_format=input_format)
    structure.to(output_file_path=output_file, output_file_format=output_format)


if __name__ == "__main__":
    convert_format()