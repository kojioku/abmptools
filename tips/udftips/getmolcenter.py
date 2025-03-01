import numpy as np


def getposmol(uobj, indexMol):
    '''
    get position of a molecule

    Args:
    uobj: universe object
    indexMol: index of molecule

    Returns:
    position: position of the molecule
    '''

    i = indexMol
    list = uobj.get("Structure.Position.mol[" + str(i) + "].atom[]")
    position = np.array(list)

    return position


def getCenter(posVec):
    '''
    get center of a molecule

    Args:
    posVec: position vector of a molecule

    Returns:
    center: center of the molecule
    '''

    center = np.average(posVec, 0)
    return center


if __name__ == "__main__":
    molid = 0

    pos = getposmol(_udf_, molid)
    cocs = getCenter(pos).tolist()

    print('Center of molecule ' + str(molid) + ':', str(cocs))
