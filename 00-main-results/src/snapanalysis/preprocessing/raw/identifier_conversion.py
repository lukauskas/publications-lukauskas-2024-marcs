# Map from the N- based id's to the new H-based IDs.
# The idea was to rename them to something that made more sense (at that point):
# assigning each histone octamer an unique H* number, and then adding M suffix
# but since then we sort of forgot why we cared about this
# nevertheless we kept the new identifiers because we got used to them.

PULLDOWN_ID_MAP = {'N11': 'H05', 'N12': 'H06', 'N13': 'H03', 'N14': 'H04', 'N15': 'H39',
                   'N16': 'H07', 'N17': 'H08', 'N18': 'H36', 'N19': 'H11', 'N2': 'H01',
                   'N20': 'H12', 'N21': 'H09', 'N22': 'H10', 'N23': 'H13', 'N24': 'H14',
                   'N25': 'H40', 'N26': 'H15', 'N27': 'H19', 'N28': 'H17', 'N29': 'H21',
                   'N30': 'H16', 'N31': 'H18', 'N32': 'H20', 'N33': 'H23', 'N34': 'H25',
                   'N35': 'H22', 'N36': 'H24', 'N37': 'H26', 'N38': 'H37', 'N39': 'H41',
                   'N40': 'H42', 'N41': 'H43', 'N42': 'H45', 'N43': 'H38', 'N44': 'H44',
                   'N45': 'H01M', 'N46': 'H03M', 'N47': 'H04M', 'N48': 'H27M', 'N49': 'H39M',
                   'N5': 'H02', 'N50': 'H07M', 'N51': 'H47', 'N52': 'H46', 'N53': 'H47M',
                   'N54': 'H46M', 'N55': 'H08M', 'N57': 'H28', 'N58': 'H29', 'N59': 'H32',
                   'N60': 'H30', 'N61': 'H31', 'N62': 'H35', 'N63': 'H34', 'N64': 'H33'}

def convert_to_new_style_identifier(old_identifier, silent=False):
    try:
        return PULLDOWN_ID_MAP[old_identifier]
    except KeyError:
        if silent:
            return None
        else:
            raise
