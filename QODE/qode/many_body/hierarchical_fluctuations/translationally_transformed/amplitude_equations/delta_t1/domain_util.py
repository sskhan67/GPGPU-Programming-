def intersection(*domains):
    if not all(domains):  return []    # relies on empty list evaluating to False when checked for truth
    else:
        max_min = max(d[0]  for d in domains)
        min_max = min(d[-1] for d in domains)
        if max_min<=min_max:  return list(range(max_min, min_max+1))
        else:                 return []

def domain(L):
    def _domain(x):  return list(range(x-L, x+L+1))
    return _domain
