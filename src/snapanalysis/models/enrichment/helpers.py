
def filter_zero_rows(matrix, min_nonzero=1):
    """
    Leaves only matrix rows hat have at least `min_nonzero` columns that are non ero
    :param matrix:
    :param min_nonzero:
    :return:
    """
    return matrix[(matrix != 0).sum(axis=1) >= min_nonzero]

def remove_all_zero_rows(matrix):
    """
    Removes all matrix rows that have all columns set to zero.

    :param matrix:
    :return:
    """
    return filter_zero_rows(matrix, min_nonzero=1)


def collapse_direction(matrix, axis=1):
    """
    Collapses matrix columns from MultiIndex (that is tricky to deal with)
    to a single-level index with F/R appended

    :param matrix:
    :param axis:
    :return:
    """

    if axis == 1:
        old_index = matrix.columns
    elif axis == 0:
        old_index = matrix.index
    else:
        raise Exception('Unsupported axis {!r}'.format(axis))

    matrix = matrix.copy()

    new_index = ['{}-{}'.format(x[1], x[0][0]) for x in old_index]

    if axis == 1:
        matrix.columns = new_index
    elif axis == 0:
        matrix.index = new_index
    else:
        raise Exception('Unsupported axis {!r}'.format(axis))

    return matrix
