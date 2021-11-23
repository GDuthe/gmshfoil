from gmsh_foil import GMSHFoil
import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input-dat-file', '-i', help="input file containing airfoil coordinates", type= str, default = None)
parser.add_argument('--foil-number-string', '-f', help="the foil number", type= str, default = '4812')
parser.add_argument('--output-mesh-file-name', '-o',help='the file name of the mesh output (without the .msh suffix)', type = str, default = None)
parser.add_argument('--view','-v',help='view the mesh after creating it',type=bool)

def _gf_mesh_run(input_dat_file, foil_number_string, output_mesh_file, view):
    g = GMSHFoil(foil_dat_file=input_dat_file, naca_foil=foil_number_string, mesh_name=output_mesh_file)
    g.create_2d_unstructured_foil_mesh()
    if view:
        g.view()


if __name__ == '__main__':
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    
    _gf_mesh_run(
            args.input_dat_file,
            args.foil_number_string,
            args.output_mesh_file_name, 
            args.view)
