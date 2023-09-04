from rdkit import Chem
from rdkit.Chem import AllChem
import os
from InchiGen import *

def GenerateSDFFromTxt(InputFile, OutputFolder, inp_type):
    Mols = []
    cwd = os.getcwd()
    InputFile
    if os.path.exists(os.path.join(cwd, InputFile)):
        fullf = os.path.join(cwd, InputFile)
    elif inp_type == 'InChI':
        fullf = os.path.join(cwd, InputFile + ".inchi")
    elif inp_type == 'Smiles':
        fullf = os.path.join(cwd, InputFile + ".smiles")
    elif inp_type == 'Smarts':
        fullf = os.path.join(cwd, InputFile + ".smarts")
    f = open(fullf, "r")
    if inp_type == 'InChI':
        for line in f.readlines():
            if line.strip() == '':
                continue
            else:
                try:
                    m = Chem.MolFromInchi(line.strip(),sanitize = True)
                    Mols.append(m)
                except:
                    print(line + " could not be read")
    elif inp_type == 'Smiles':
        for line in f.readlines():
            if line.strip() == '':
                continue
            else:
                try:
                    m = Chem.MolFromSmiles(line.strip(),sanitize = True)
                    Mols.append(m)
                except:
                    print(line + " could not be read")
    elif inp_type == 'Smarts':
        for line in f.readlines():
            if line.strip() == '':
                continue
            else:
                try:
                    m = Chem.MolFromSmarts(line.strip(),sanitize = True)
                    Mols.append(m)
                except:
                    print(line + " could not be read")
    else:
        print('unrecognised')
    GeneratedFiles = GenerateSDFFromMols(Mols, OutputFolder, inp_type)
    f.close()
    return GeneratedFiles

#===========================

def GenerateSDFFromMols(mols, OutputFolder, inp_type):
    files = []
    cwd = os.getcwd()
    for ind, m in enumerate(mols):
        m = AllChem.AddHs(m, addCoords=True)
        AllChem.EmbedMolecule(m)
        AllChem.MMFFOptimizeMolecule(m)
        Chem.rdmolops.AssignStereochemistryFrom3D(m)
        f = inp_type +  '_Mol_' + str(ind) +"_.sdf"
#        files.append(f[:-4])
#        fullf = os.path.join(cwd, f)
        fullf = os.path.join(OutputFolder, f)
#        sdfFile = open(fullf, 'a')
        sdfFile = open(fullf, 'w')
        save3d = Chem.SDWriter(sdfFile)
        save3d.write(m)
        files.append(fullf)
    return files

#===========================

def GenerateMolFromSDF(InputFile):
    cwd = os.getcwd()
    fullf = os.path.join(cwd, InputFile + '.sdf')
    m = Chem.MolFromMolFile(fullf, removeHs=False,sanitize = True)
    return m

#===========================

def CleanUp(InputFiles, OutputFolder):
    # check input file types
    CleanedInputFiles = []
    cwd = os.getcwd()
#    fullf = os.path.join(OutputFolder, f + 'cleaned.sdf')
    fullf = os.path.join(OutputFolder, 'Input_cleaned.sdf')
    for f in InputFiles:
        if f.endswith('.sdf'):
            f = f[:-4]
        m = GenerateMolFromSDF(f)
        m = AllChem.AddHs(m, addCoords=True)
        AllChem.EmbedMolecule(m)
        AllChem.MMFFOptimizeMolecule(m)
        Chem.rdmolops.AssignStereochemistryFrom3D(m)
        sdfFile = open(fullf, 'a')
        save3d = Chem.SDWriter(sdfFile)
        save3d.write(m)

    CleanedInputFiles.append(fullf)
    return CleanedInputFiles

#===========================

def NumberofStereoCentres(InputFile):
    cwd = os.getcwd()
#    fullf = os.path.join(cwd, InputFile + '.sdf')
    fullf = os.path.join(cwd, InputFile)
    m = Chem.MolFromMolFile(fullf, removeHs=False)
    m = AllChem.AddHs(m, addCoords=True)
    AllChem.EmbedMolecule(m)
    AllChem.MMFFOptimizeMolecule(m)
    Chem.rdmolops.AssignStereochemistryFrom3D(m)
    nStereo = Chem.rdMolDescriptors.CalcNumAtomStereoCenters(m)
#    print (m)
    return nStereo

#===========================


