from family import T_family, transit_M_dedup
bad = [b for b in range(6, 101, 2) if transit_M_dedup(T_family(b)) != 3]
print('even b in [6,100] with M != 3:', bad if bad else 'NONE - all equal 3')
