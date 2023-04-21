#!/Users/mac/opt/anaconda3/envs/research/bin/python3
"""
Usage
-----
    $ python3 addMag.py  --atom_config /Users/mac/我的文件/Mycode/new/new2/matersdk/test_data/atom_config/addMag/output/atom.config --atom_config_template /Users/mac/我的文件/Mycode/new/new2/matersdk/test_data/atom_config/addMag/tem/atom.config -n 12 -n 25
"""
import click
from masksdk.io.pwmat.addMag import AtomConfigMagTem


@click.command()
@click.option('--atom_config',
            #default=1,
            type=str,
            prompt='The path of atom.config',
            help="he path of atom.config.",
            )
@click.option('--atom_config_template',
            #default=1,
            type=str,
            prompt='The path of atom.config (template)',
            help="he path of atom.config (template).",
            )
@click.option('--atomic_number_zero_lst', '-n',
            #default=[12],
            prompt='The atomic numbers whose magnetic moment is 0',
            help="The atomic numbers whose magnetic moment is 0.",
            multiple=True,
            )
def add_mag(atom_config, atom_config_template, atomic_number_zero_lst):
    atom_config_mag_tem = AtomConfigMagTem(
                            atom_config_path_sqs=atom_config,
                            atom_config_path_template=atom_config_template
                            )
    atomic_number_zero_lst = [int(value) for value in atomic_number_zero_lst]
    atom_config_mag_tem.assign_magnetic_moment(
                            atomic_number_zero_lst=atomic_number_zero_lst,
                            )
    
    atom_config_mag_tem.append_mag_to_atomconfig()


if __name__ == "__main__":
    add_mag()