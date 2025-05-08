def kmph_2_mps(v_kmph):
    """
    Convert speed from kilometers per hour to meters per second.
    
    Parameters:
    v_kmph (float): Speed in kilometers per hour.
    
    Returns:
    float: Speed in meters per second.
    """
    return v_kmph / 3.6

def kg_2_N(kg):
    """
    Convert mass from kilograms to newtons.
    
    Parameters:
    kg (float): Mass in kilograms.
    
    Returns:
    float: Mass in newtons.
    """
    return kg * 9.81  # 1 kg = 9.81 N (acceleration due to gravity)

def N_2_kg(N):
    """
    Convert force from newtons to kilograms.
    
    Parameters:
    N (float): Force in newtons.
    
    Returns:
    float: Force in kilograms.
    """
    return N / 9.81  # 1 N = 1 kg * 9.81 m/s^2 (acceleration due to gravity)