from geopy.geocoders import Nominatim

geolocator = Nominatim(user_agent=f"{pays}continentcharles.seanna@outlook.fr")
location = geolocator.geocode(pays)


latitude=location.latitude
longitude=location.longitude
if -35 <= latitude <= 37 and -20 <= longitude <= 55:
    return( "Africa")
elif 35 <= latitude <= 70 and -10 <= longitude <= 40:
    return( "Europe")
elif 5 <= latitude <= 55 and 70 <= longitude <= 140:
    return( "Asia")
elif -55 <= latitude <= 15 and -80 <= longitude <= -35:
    return( "South America")
elif 15 <= latitude <= 72 and -170 <= longitude <= -50:
    return( "North America")
elif -50 <= latitude <= -10 and 110 <= longitude <= 180:
    return( "Oceania")
else:
    return(pays)