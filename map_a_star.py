import csv
import folium
import pandas as pd
import os
import webbrowser
from a_star import *

def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ') 
    directory_read = os.path.join("Outputs" ,  date) 
    return date,directory_read

def readcitiesOutput(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    # Create a map centered on the mean latitude and longitude
    m = folium.Map(location=[df['latitude'].mean(), df['longitude'].mean()], zoom_start=2)
    
    start_city = int(input('Starting from: '))
    end_city = int(input('Going to: '))
    readA_starFuelPriceOutput(date, directory_read)
    readA_starvelocityOutput(date, directory_read)
    readA_starfuelConsumptionOutput(date, directory_read)
    readA_starTollPriceOutput(date, directory_read)
    readA_star_distancesOutput(date, directory_read)
    readA_starTimeOutput(date, directory_read)
    readA_startotalCostOutput(date, directory_read)
    
    path, A_stardistances, A_starTime, A_starTollPrice, A_startotalCost = sol_cost_with_ids(start_city,end_city)
    for _, row in df.iterrows():
        color = 'blue'  # set the default marker color to blue
        if row['City ID'] == start_city:
            color = 'green'  # set the marker color to green for the start city
        elif row['City ID'] == end_city:
            color = 'red'  # set the marker color to red for the end city
        elif row['City ID'] in path:
            color = 'purple'  # set the marker color to purple for cities in the path


        folium.Marker(
            location=[row['latitude'], row['longitude']],
            popup=int(row['City ID']),
            icon=folium.Icon(color=color)
        ).add_to(m)

