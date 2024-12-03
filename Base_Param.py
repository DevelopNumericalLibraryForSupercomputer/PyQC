def unit_conversion(string,unit1,unit2):
    eng_au2ev = float(27.211324570273)

    if string=='energy':
       if unit1=='au' and unit2=='ev':
          return eng_au2ev
       elif unit1=='ev' and unit2=='au':
          return 1.0/eng_au2ev 

    str0 = str(string)+' : '+str(unit1)+' -> '+str(unit2)
    print(" * Error: No string found : "+str0)
    return 0.0

