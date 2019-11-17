import ephem
import numpy as np

def localtime_to_utc(local,mdt=True):
    """

    Converts local time at APO to UTC

    Parameters
    ----------
    local : str
        Local time and date given in the following format: 'YYYY/MM/DD H:M:S.'
    mdt : bool, optional
        Mountain Daylight time. Default is True

    Returns
    -------
    utc : str
        Local time and date converted to Coordinated Universal Time.
    """

    #Separate time and date
    time = local.splt(' ')[0]
    date = local.split(' ')[1]

    #Separate hour and minute
    hour = int(time.split(':')[0])
    minute = int(time.split(':')[1])
    second = int(time.split(':')[-1])

    #separate the day
    year = int(date.split('/')[0])
    month = int(date.split('/')[1])
    day = int(date.split('/')[-1])

    #Daylight savings time makes it UTC-6:00
    if mdt:
        hour = hour + 6

        if hour > 24:
            day = day + 1
            hour = hour%24
    else:
        hour = hour + 7

        if hour > 24:
            day = int(day) + 1
            hour = hour%24

    return "{}/{}/{} {}:{}:{}".format(year,month,day,hour,minute,second)

def localtime_to_lst(local,
                     mdt=True,
                     latitude='32.7803',
                     longitude='-108.8203'
                     ):
    """

    Converts local time at APO to local sidereal time

    Parameters
    ----------
    local : str
        Local time and date given in the following format: 'YYYY/MM/DD H:M:S.'
    mdt : bool, optional
        Mountain Daylight time. Default is True
    latitude : str
        Latitude of site. Set to 32.7803
    longitude : str
        Longitude of site. Set to -108.8203

    Returns
    -------
    lst : str
        Local time converted to local sidereal time.
    """

    #Separate time and date
    time = local.splt(' ')[0]
    date = local.split(' ')[1]

    #Separate hour and minute
    hour = int(time.split(':')[0])
    minute = int(time.split(':')[1])
    second = int(time.split(':')[-1])

    #separate the day
    year = int(date.split('/')[0])
    month = int(date.split('/')[1])
    day = int(date.split('/')[-1])

    #Daylight savings time makes it UTC-6:00
    if mdt:
        hour = hour + 6

        if hour > 24:
            day = day + 1
            hour = hour%24
    else:
        hour = hour + 7

        if hour > 24:
            day = int(day) + 1
            hour = hour % 24

    utc = "{}/{}/{} {}:{}:{}".format(year,month,day,hour,minute,second)

    #ephem is a great package for getting things like local sidereal time, coordinate conversions, etc.
    apo = ephem.Observer()
    apo.lat, apo.lon = latitude, longitude
    apo.date = utc

    #print the sidereal_time
    print(apo.sidereal_time())


def equatorial_to_horizontal(declination,
                             right_ascension,
                             local_sidereal_time,
                             observers_latitude = '32.7803'
                             ):


    """

    Converts from equatorial coordinate system to horizontal coordinate system.
    
    Parameters
    ----------
    declination : str
        Declination of the target in sexagesimal format 'dd:mm:ss.s'.
    right_ascension : str
        Right ascension of the target in sexagesimal format 'hh:mm:ss'.
    local_sidereal_time : str
        The local sidereal time which you wish to compute the horizontal coordinates for. Given in Military format.
    observers_latitude : str, optional
        The altitude of the observatory. Defaults to 32.7803.

    Return
    -------
    az, alt : tuple
        The azimuth and altitude based off the equatorial coordinates provided.
    """

    #Split the Right Ascension into it's components.
    ra_hour = right_ascension.split(':')[0]
    ra_minute = right_ascension.split(':')[1]
    ra_second = right_ascension.split(':')[2]

    #Convert Right Ascension into radians.
    ra = (int(ra_hour) + int(ra_minute)/60. + int(ra_second)/3600.) * (np.pi/12.)


    #Repeat for declination
    dec_hour = declination.split(':')[0]
    dec_minute = declination.split(':')[1]
    dec_second = declination.split(':')[-1]

    dec = (int(dec_hour)+int(dec_minute)/60.+int(dec_second)/3600.) *(np.pi/180)

    #Repeat for local sidereal time
    lst_hour = local_sidereal_time.split(':')[0]
    lst_minute = local_sidereal_time.split(':')[1]
    lst_second = local_sidereal_time.split(':')[2]

    lst = (int(lst_hour)+int(lst_minute)/60.+int(lst_second)/3600.)*(np.pi/12)

    #And lastly for the latitude
    lat = observers_latitude*(np.pi/180.)

    #Calculate altitude (Eq. 1.76)
    sin_alt = np.sin(lat) * np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(lst - ra)

    alt = np.asin(sin_alt)
    alt_deg = np.degrees(alt)

    #Calculate azimuth Eq. (1.77)
    sin_az = np.cos(dec)*np.sin(lst - ra)/np.cos(alt)

    # Eq. (1.78)
    cos_az = (np.sin(lat) * np.cos(dec) * np.cos(lst - ra) - np.cos(lat) * np.sin(dec)) / np.cos(alt)

    #Take care of quadrant ambiguity.
    az = np.arctan2(sin_az, cos_az)
    az_deg = mnpdegrees(az) % 360

    return az_deg, alt_deg


def dec_to_sexa(decimal, key, round_value=1):
    '''
    Inputs:
        Decimal is the floating value of the number you want to convert
        into sexagesimal form. It must either be in hours or degrees.
        Key is a string, either "h" for hour or "d" for degrees.
        Round_value is optional, default set to 1 decimal place.

    Returns:
        The number in sexagesimal form
    '''
    dec_unit = int(decimal)
    dec_minute = (abs(decimal) - abs(int(decimal))) * 60
    dec_sec = (dec_minute - int(dec_minute)) * 60

    return "{}{} {}m {}s".format(dec_unit, key, int(dec_minute), round(dec_sec, round_value))