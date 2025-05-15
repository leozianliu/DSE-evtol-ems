import csv
import requests
import time

# Function to geocode an address using Nominatim (OpenStreetMap)
def geocode_address(hospital_name):
    # Add "Bavaria, Germany" to improve search accuracy
    search_term = f"{hospital_name}, Bayern, Deutschland"
    
    # OpenStreetMap Nominatim API with appropriate user-agent
    url = "https://nominatim.openstreetmap.org/search"
    headers = {"User-Agent": "HospitalLocationFinder/1.0"}
    params = {
        "q": search_term,
        "format": "json",
        "limit": 1
    }
    
    try:
        response = requests.get(url, headers=headers, params=params)
        data = response.json()
        
        if data and len(data) > 0:
            lat = float(data[0]["lat"])
            lon = float(data[0]["lon"])
            return lat, lon
        else:
            return None, None
    except Exception as e:
        print(f"Error geocoding {hospital_name}: {e}")
        return None, None
    finally:
        # Be respectful of the API usage policy
        time.sleep(1)

# List of hospitals that likely operate 24/7 (major medical centers, university hospitals)
# This is a pre-filtered list focusing on major facilities likely to be 24-hour
likely_24hr_hospitals = [
    "Universitätsklinikum Augsburg - Medizincampus",
    "Universitätsklinikum Würzburg",
    "Universitätsklinikum Erlangen",
    "LMU Klinikum München - Campus Großhadern",
    "Klinikum Nürnberg Nord",
    "Klinikum rechts der Isar der Technischen Universität München - Stammgelände",
    "Klinikum Ingolstadt",
    "München Klinik Bogenhausen, München Klinik gGmbH",
    "Krankenhaus Barmherzige Brüder Regensburg",
    "Klinikum Nürnberg Süd",
    "Universitätsklinikum Regensburg",
    "Klinikum Fürth",
    "Klinikum Aschaffenburg-Alzenau gemeinnützige GmbH",
    "München Klinik Harlaching, München Klinik gGmbH",
    "Klinikum Bayreuth",
    "Klinikum Passau",
    "RoMed Klinikum Rosenheim",
    "BG Klinikum Murnau gGmbH",
    "Klinikum Traunstein",
    "München Klinik Neuperlach, München Klinik gGmbH",
    "Klinikum St. Marien Amberg",
    "Klinikum Landshut AdöR der Stadt Landshut",
    "Klinikverbund Allgäu gGmbH - Klinikum Kempten",
    "Klinikum Memmingen AöR",
    "DONAUISAR Klinikum Deggendorf",
    "Klinikum St. Elisabeth Straubing GmbH",
    "Sana Klinikum Hof",
    "Helios Amper-Klinikum Dachau",
    "Rotkreuzklinikum München gGmbH",
    "Helios Klinikum München West",
    "InnKlinikum Altötting",
    "Barmherzige Brüder Krankenhaus München",
    "Klinikum Freising GmbH",
    "Krankenhaus Agatharied",
    "Caritas-Krankenhaus St. Josef",
    "Krankenhaus Martha-Maria Nürnberg",
    "Klinikum Starnberg",
    "Malteser Waldkrankenhaus St. Marien",
    "Klinikum Landkreis Erding",
    "Helios Frankenwaldklinik Kronach",
    "InnKlinikum Altötting und Mühldorf / InnKlinikum Mühldorf",
    "Krankenhaus Rummelsberg",
    "Kreisklinik Roth",
    "Krankenhaus St. Barbara (Schwandorf)",
    "Donau-Ries Klinik Donauwörth",
    "Sana Kliniken des Landkreises Cham - Krankenhaus Cham",
    "Klinikum Landsberg am Lech",
    "Klinikum Altmühlfranken Gunzenhausen",
    "Deutsches Herzzentrum München des Freistaates Bayern"
]

def main():
    # Output file
    with open('bavaria_24hr_hospitals_coordinates.csv', 'w', newline='', encoding='utf-8') as outfile:
        writer = csv.writer(outfile)
        writer.writerow(['Hospital Name', 'Latitude', 'Longitude'])
        
        for hospital in likely_24hr_hospitals:
            print(f"Geocoding: {hospital}")
            lat, lon = geocode_address(hospital)
            
            if lat and lon:
                writer.writerow([hospital, lat, lon])
                print(f"Found coordinates: {lat}, {lon}")
            else:
                writer.writerow([hospital, "Not found", "Not found"])
                print("Coordinates not found")
            
            # Small delay to be nice to the API
            time.sleep(1)

if __name__ == "__main__":
    main()