from random_instance import *
import ast
import math
import os
import pandas


def wanted_date():
    options = os.listdir(os.path.join("Outputs", ))
    print(options)

    date = input('date: ')
    directory_read = os.path.join("Outputs", date)
    return date, directory_read

A_stardistances = []
A_starcities = []
A_startotalCost = []
A_starFuelConsumption = []
A_starFuelPrice = []
A_starTime = []
A_starTollPrice =[]
A_starVelocity = [] 
fuel_needed = []



def readA_star_distancesOutput(date, directory_read):
    distances_file_name = 'distances_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, distances_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    distances_read = df['distances'].tolist()

    for i in range(len(distances_read)):
        A_stardistances.append([from_read[i], to_read[i], distances_read[i]])

def readA_star_citiesOutput(date, directory_read):
    cities_file_name = 'cities_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, cities_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    cityId_read = df['City ID'].tolist()
    x_read = df['latitude'].tolist()
    # extract the 'time' column as a list
    y_read = df['longitude'].tolist()

    for i in range(len(cityId_read)):
        A_starcities.append([cityId_read[i], x_read[i], y_read[i]])

def readA_startotalCostOutput(date, directory_read):
    totalCost_file_name = 'totalCost_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, totalCost_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'totalCost' column as a list
    totalCost_read = df['Total Cost'].tolist()

    for i in range(len(totalCost_read)):
        A_startotalCost.append([from_read[i], to_read[i], totalCost_read[i]])

def readA_starfuelConsumptionOutput(date, directory_read):
    FuelConsumption_file_name = 'fuel_consumption_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, FuelConsumption_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    FuelConsumption_read = df['FuelConsumption'].tolist()

    for i in range(len(FuelConsumption_read)):
        A_starFuelConsumption.append([from_read[i], to_read[i], FuelConsumption_read[i]])

def readA_starFuelPriceOutput(date, directory_read):
    fuelPrice_file_name = 'fuel_price_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, fuelPrice_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    City_ID_read = df['cityID'].tolist()
    # extract the 'time' column as a list
    fuelPrice_read = df['fuelPrice'].tolist()

    for i in range(len(fuelPrice_read)):
        A_starFuelPrice.append([City_ID_read[i], fuelPrice_read[i]])

def readA_starTimeOutput(date, directory_read):
    time_file_name = 'time_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, time_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    time_read = df['time'].tolist()

    for i in range(len(time_read)):
        A_starTime.append([from_read[i], to_read[i], time_read[i]])

def readA_starTollPriceOutput(date, directory_read):
    TollPrice_file_name = 'toll_price_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, TollPrice_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    TollPrice_read = df['TollPrice'].tolist()

    for i in range(len(TollPrice_read)):
        A_starTollPrice.append([from_read[i], to_read[i], TollPrice_read[i]])

def readA_starvelocityOutput(date, directory_read):
    velocity_file_name = 'velocity_' + date + '.csv'

    # construct the full file path using the 'directory' variable
    file_path = os.path.join(directory_read, velocity_file_name)

    # read the CSV file into a DataFrame
    df = pd.read_csv(file_path, sep=';')
    from_read = df['from'].tolist()
    to_read = df['to'].tolist()
    # extract the 'time' column as a list
    velocity_read = df['velocity'].tolist()

    for i in range(len(velocity_read)):
        A_starVelocity.append([(from_read[i], to_read[i]), velocity_read[i]])

def finding_fuel_needed():
    for i in range(len(A_stardistances)):
        fuel_needed.append([A_stardistances[i][0] , A_stardistances[i][1] , A_stardistances[i][2] * A_starFuelConsumption[i][2]])

#call the readAllAsCSV function before using the a_star function to populate the required data from the CSV files.
def readAllAsCSV():
    date, directory_read = wanted_date()
    readA_star_distancesOutput(date, directory_read)
    readA_star_citiesOutput(date, directory_read)
    readA_startotalCostOutput(date, directory_read)
    readA_starfuelConsumptionOutput(date, directory_read)
    readA_starFuelPriceOutput(date, directory_read)
    readA_starTimeOutput(date, directory_read)
    readA_starTollPriceOutput(date, directory_read)
    readA_starvelocityOutput(date, directory_read)
    finding_fuel_needed()


#Adding heuristic to find the estimated cost of travelling from the current city to the end city (cumulated), heuristic needs to be added before working on a_star search

#Heuristic distance: Manhattan or Euclidian distance applications will be searched for min weight distance heuristic

def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the distance between two points on the Earth's surface using the Haversine formula.
    """
    # Radius of the Earth in kilometers
    radius = 6371

    # Convert latitude and longitude from degrees to radians
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    # Calculate the differences between the latitudes and longitudes
    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad

    # Calculate the square of half the chord length between the points
    a = math.sin(delta_lat / 2) ** 2 + math.cos(lat1_rad) * math.cos(lat2_rad) * math.sin(delta_lon / 2) ** 2

    # Calculate the angular distance in radians
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    # Calculate the distance
    distance = radius * c
    return distance

def manhattan(lat1, lon1, lat2, lon2):
    """
    Calculate the Manhattan distance between two points on the Earth's surface.
    """
    # Convert latitude and longitude from degrees to kilometers
    lat_km = (lat2 - lat1) * 111
    lon_km = (lon2 - lon1) * 111

    # Calculate the Manhattan distance
    distance = abs(lat_km) + abs(lon_km)
    return distance

def euclidean(lat1, lon1, lat2, lon2):
    """
    Calculate the Euclidean distance between two points on the Earth's surface.
    """
    # Convert latitude and longitude from degrees to kilometers
    lat_km = (lat2 - lat1) * 111
    lon_km = (lon2 - lon1) * 111

    # Calculate the Euclidean distance
    distance = math.sqrt(lat_km ** 2 + lon_km ** 2)
    return distance

def get_city_by_id(city_id, cities):
    """
    Get the city data by its ID.
    """
    for city in cities:
        if city['id'] == city_id:
            return city
    return None

def heuristic_distance_haversine(city_id, end_city_id, cities):
    current_city = get_city_by_id(city_id, cities)
    end_city = get_city_by_id(end_city_id, cities)
    distance = haversine(
        current_city['latitude'],
        current_city['longitude'],
        end_city['latitude'],
        end_city['longitude']
    )
    return distance

def heuristic_distance_manhattan(city_id, end_city_id, cities):
    current_city = get_city_by_id(city_id, cities)
    end_city = get_city_by_id(end_city_id, cities)
    distance = manhattan(
        current_city['latitude'],
        current_city['longitude'],
        end_city['latitude'],
        end_city['longitude']
    )
    return distance

def heuristic_distance_euclidean(city_id, end_city_id, cities):
    current_city = get_city_by_id(city_id, cities)
    end_city = get_city_by_id(end_city_id, cities)
    distance = euclidean(
        current_city['latitude'],
        current_city['longitude'],
        end_city['latitude'],
        end_city['longitude']
    )
    return distance




#A-star algorithm steps to follow:
#1. Set all point distance to infinity except for the starting point set distance to 0. +
#2. Set all points, including start point as a non-visited node.  +
#3. Set the non-visited node with the smallest current distance as the current node "C." 
#4. For each neighbor “N” of your current node: add the current distance of “C” with the weight of the edge connecting “C” – “N” and the weight to the destination point (heuristic). If it's smaller than the current distance of “N," set it as the new current distance of “N."
#5. Mark the current node “C” as visited. 
#6. Repeat the step above from step 3 until one of the neighbors “N” is the destination point.
#Following those steps one by one:


#1. Set all point distance to infinity except for the starting point set distance to 0. The distance to the start city is set to 0 in the min_costs dictionary.


#This function takes into account two main things which are distances and the costs:
def A_star_path(start_city_id, end_city_id, costs, distances, cities, distance_method):
 
 
 #Initially, all cities except the start city have a distance of infinity (not explicitly set in the code) since they haven't been visited yet. 
 min_costs = {start_city_id: 0}
 
 #Dictionary created to store the previous city in the optimal path.
 previous_city = {}
 
 #Created a set to keep track of visited cities.
 visited = set()
 
 # Start with the start_city as the current city
 current_city = start_city_id 

 if distance_method == 'euclidean':
     distance_function = heuristic_distance_euclidean
 if distance_method == 'manhattan':
     distance_function = heuristic_distance_manhattan
 if distance_method == 'haversine':
     distance_function = heuristic_distance_haversine
 else:
     raise ValueError("Invalid distance calculation method. Please select between euclidian, manhattan or haversine")

    
 while current_city != end_city_id:
        # Mark the current city as visited
        visited.add(current_city)

        # Get the cost to reach the current city
        current_cost = min_costs[current_city]

        # Find the neighbors of the current city, 
        neighbors = [c2 for (c1, c2, _) in costs if c1 == current_city]
        
        for neighbor in neighbors:
            # Check if the distance is -1, skip this city
            if any(c1 == current_city and c2 == neighbor and d == -1 for (c1, c2, d) in distances):
                continue
            # Calculate the distance to reach the neighbor from the start_city
            neighbor_distance = current_cost + next(c for (c1, c2, c) in distances if c1 == current_city and c2 == neighbor)
            
            #Checks if the neighbor city is not already in the min_costs dictionary or if the newly calculated cost (neighbor_cost) is smaller than the previously recorded cost for the neighbor
            if neighbor not in min_costs or neighbor_distance < min_costs[neighbor]:
                # Update the minimum cost and previous city for the neighbor
                min_costs[neighbor] = neighbor_distance
                previous_city[neighbor] = current_city
        # Find the next unvisited city with the minimum cost
        unvisited_cities = {city: min_costs[city] for city in min_costs if city not in visited}

        if not unvisited_cities:
            # No unvisited cities remaining, end_city_id is not reachable
            raise ValueError("End city is not reachable from the start city.")
        
        # Calculate the heuristic value for each unvisited city (added options of haversine, euclidian and manhattan)
        heuristic_values = {city: distance_function(city, end_city_id, cities) for city in unvisited_cities}

        
        current_city = min(heuristic_values, key=heuristic_values.get)

 #Recounstruct the optimal path
 path = []
 city = end_city_id
 while city != start_city_id:
     path.append(city)
     city = previous_city[city]
 path.append(start_city_id)

 #Reverse the path for building it backwards    
 path.reverse()
 return path

def A_star_distance(start_city_id, end_city_id, distances, heuristic_function):
    # Initially, all cities except the start city have a distance of infinity
    min_distances = {start_city_id: 0}
    
    # Dictionary to store the previous city in the optimal path
    previous_city = {}
    
    # Set to keep track of visited cities
    visited = set()
    
    # Start with the start_city as the current city
    current_city = start_city_id
    
    while current_city != end_city_id:
        # Mark the current city as visited
        visited.add(current_city)
        
        # Get the distance to reach the current city
        current_distance = min_distances[current_city]
        
        # Find the neighbors of the current city
        neighbors = [c2 for (c1, c2, _) in distances if c1 == current_city]
        
        for neighbor in neighbors:
            # Check if the distance is -1, skip this city
            if any(c1 == current_city and c2 == neighbor and d == -1 for (c1, c2, d) in distances):
                continue
            
            # Calculate the distance to reach the neighbor from the start_city
            neighbor_distance = current_distance + next(c for (c1, c2, c) in distances if c1 == current_city and c2 == neighbor)
            
            if neighbor not in min_distances or neighbor_distance < min_distances[neighbor]:
                # Update the minimum distance and previous city for the neighbor
                min_distances[neighbor] = neighbor_distance
                previous_city[neighbor] = current_city
        
        # Find the next unvisited city with the minimum distance plus heuristic value
        unvisited_cities = {city: min_distances[city] + heuristic_function(city, end_city_id) for city in min_distances if city not in visited}
        
        if not unvisited_cities:
            # No unvisited cities remaining, end_city_id is not reachable
            raise ValueError("End city is not reachable from the start city.")
        
        current_city = min(unvisited_cities, key=unvisited_cities.get)
    
    # Reconstruct the optimal path
    path = []
    city = end_city_id
    while city != start_city_id:
        path.append(city)
        city = previous_city[city]
    path.append(start_city_id)
    
    # Reverse the path since we built it backwards
    path.reverse()
    
    return path


def minimize_fuel_cost(path, fuel_needed, fuel_price, max_fuel_capacity, starting_fuel):
    fuel_cost = 0  # Initialize total fuel cost
    current_fuel = starting_fuel  # Initialize current fuel level
    fuel_purchases = []  # List to store fuel purchases
    
    fuel_prices_in_path = {}  # Dictionary to store fuel prices of cities in the path
    
    # Find the fuel prices of cities in the path
    for city, price in fuel_price:
        if city in path:
            fuel_prices_in_path[city] = price
    
    # Create a copy of fuel_prices_in_path
    modified_fuel_prices = fuel_prices_in_path.copy()

    # Set the last element's value to infinity
    last_key = path[-1]
    modified_fuel_prices[last_key] = float('inf')

    # Sort the modified dictionary by values
    sorted_fuel_prices = sorted(modified_fuel_prices.items(), key=lambda x: x[1])
    

       
    # Iterate through each city in the path
    for i in range(len(path)-1):
        source_city = path[i]
        destination_city = path[i+1]
        fuel_needed_current_leg = None
        
        # Find the fuel needed for the current leg of the journey
        for fuel_leg in fuel_needed:
            if fuel_leg[0] == source_city and fuel_leg[1] == destination_city:
                fuel_needed_current_leg = fuel_leg[2]
                break
        
        # Check if the destination is the last city
        if i ==  path[len(path) - 2] and i == path[0]:
            
            # Finish the journey at the destination city (no need to buy fuel)
            if fuel_needed_current_leg <= current_fuel:
                current_fuel = current_fuel - fuel_needed_current_leg
                continue
            else:
                # Fill the tank with the amount of fuel needed for the current leg
                fuel_to_buy = fuel_needed_current_leg - current_fuel
                fuel_cost += fuel_to_buy * fuel_prices_in_path[source_city]
                current_fuel = 0
                fuel_purchases.append((source_city, fuel_to_buy))
                break
                
                
                
                
        # Check if the destination is the minimum fuel price city
        if source_city == sorted_fuel_prices[0][0]:
            # Fill the tank with the amount of fuel on hand to reach the next city
            fuel_to_buy = max_fuel_capacity - current_fuel
            fuel_cost += fuel_to_buy * fuel_prices_in_path[source_city]
            current_fuel = max_fuel_capacity - fuel_needed_current_leg
            fuel_purchases.append((source_city, fuel_to_buy))
            
        else:
            # Check if it's possible to reach the destination with current fuel
            if fuel_needed_current_leg <= current_fuel:
                current_fuel = current_fuel - fuel_needed_current_leg
                # Find the next minimum fuel price city that is reachable with current fuel
                continue       
            else:
                # Fill the tank with the amount of fuel needed for the current leg
                fuel_to_buy = fuel_needed_current_leg - current_fuel
                fuel_cost += fuel_to_buy * fuel_prices_in_path[source_city]
                current_fuel = 0
                fuel_purchases.append((source_city, fuel_to_buy))
    
    # Return the list of fuel purchases and the total fuel cost
    return fuel_purchases, fuel_cost

#This function calculates the total distance travelled via looking up the distances between two cities and addind them up cumulatively until the destination city is reached.
def calculate_path_distance(path, distances):
    total_distance = 0
    
    for i in range(len(path)-1):
        city_id1 = path[i]
        city_id2 = path[i+1]
        
        # Look up the distance between the two cities
        for (city_1,city_2, distance) in distances:
            if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                total_distance += distance
                break
                
    return total_distance


def calculate_path_time(path, times):
  
    total_time = 0
    
    for i in range(len(path)-1):
        city_id1 = path[i]
        city_id2 = path[i+1]
        
        # Look up the time between the two cities
        for (city_1,city_2, time) in times:
            if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                total_time += time
                break
                
    return total_time

def calculate_path_toll(path, tolls):
    
    total_toll = 0
    
    for i in range(len(path)-1):
        city_id1 = path[i]
        city_id2 = path[i+1]
        
        # Look up the toll between the two cities
        for (city_1,city_2, toll) in tolls:
            if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                total_toll += toll
                break
                
    return total_toll

def calculate_path_total_cost(path, costs):
    
        total_total_cost = 0
        
        for i in range(len(path)-1):
            city_id1 = path[i]
            city_id2 = path[i+1]
            
            # Look up the total_cost between the two cities
            for (city_1,city_2, total_costs) in costs:
                if (city_id1, city_id2) == (city_1,city_2) or (city_id2, city_id1) == (city_1,city_2):
                    total_total_cost += total_costs
                    break
                    
        return total_total_cost
    
def sol_distances_with_ids(i, j):
    start_city_id = int(i)
    end_city_id = int(j)
    path = A_star_distance(start_city_id, end_city_id, A_stardistances)
    total_distance = calculate_path_distance(path, A_stardistances) 
    total_time = calculate_path_time(path, A_starTime) 
    total_toll = calculate_path_toll(path, A_starTollPrice)
    total_cost = calculate_path_total_cost(path, A_startotalCost)
    print(f"Path from {start_city_id} to {end_city_id}: {path}")
    print(f"Total Distance: {total_distance}")
    print(f'Total Time:  {total_time}')
    print(f'Total Toll Price: {total_toll}')
    print(f'Total Cost: {total_cost}')

    return path, total_distance, total_time, total_toll, total_cost

def sol_for_distance():
    IDs = []
    for id in range(len(A_starcities)):
        IDs.append(A_starcities[id][0])
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = A_star_distance(start_city_id, end_city_id, A_stardistances)
                total_distance = calculate_path_distance(path, A_stardistances) 
                total_time = calculate_path_time(path, A_starTime) 
                total_toll = calculate_path_toll(path, A_starTollPrice)
                total_cost = calculate_path_total_cost(path, A_startotalCost)
                print(f"Path from {start_city_id} to {end_city_id}: {path}")
                print(f"Total Distance: {total_distance}")
                print(f'Total Time:  {total_time}')
                print(f'Total Toll Price: {total_toll}')
                print(f'Total Cost: {total_cost}')

allA_star_distance_output = []
def save_for_distance():

    IDs = []
    for id in range(len(A_starcities)):
        IDs.append(A_starcities[id][0])
    
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = A_star_distance(start_city_id, end_city_id, A_stardistances)
                total_distance = calculate_path_distance(path, A_stardistances) 
                total_time = calculate_path_time(path, A_starTime) 
                total_toll = calculate_path_toll(path, A_starTollPrice)
                total_cost = calculate_path_total_cost(path, A_startotalCost)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, A_starFuelPrice, F, Sf)
                allA_star_distance_output.append([start_city_id, end_city_id, total_cost + fuel_cost, total_distance, path, total_toll, total_time, fuel_purchase, fuel_cost])
                
def saveDistanceSolutionOutput():
    df_SolutionDistMatrix = pd.DataFrame(allA_star_distance_output)
    # define csv file name
    solD_file_name = 'solution_distance_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_SolutionDistMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City_Amount', 'fuel cost']

    if not os.path.exists(sol_directory):
        os.makedirs(sol_directory)
    # writing data frame to a CSV file
    df_SolutionDistMatrix.to_csv(os.path.join(sol_directory, solD_file_name), sep=';', index=False)
   
def sol_cost_with_ids(i,j):
    
    start_city_id = int(i)
    end_city_id = int(j)
    path = A_star_path(start_city_id, end_city_id, A_startotalCost, A_stardistances)
    total_distance = calculate_path_distance(path, A_stardistances) 
    total_time = calculate_path_time(path, A_starTime) 
    total_toll = calculate_path_toll(path, A_starTollPrice)
    total_cost = calculate_path_total_cost(path, A_startotalCost)
    
    return path, total_distance, total_time, total_toll, total_cost

def sol_for_cost():
    IDs = []
    for id in range(len(A_starcities)):
        IDs.append(A_starcities[id][0])
    
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = A_star_path(start_city_id, end_city_id, A_startotalCost, A_stardistances)
                total_distance = calculate_path_distance(path, A_stardistances) 
                total_time = calculate_path_time(path, A_starTime) 
                total_cost = calculate_path_total_cost(path, A_startotalCost)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, A_starFuelPrice, F, Sf)
                print(f"Path from {start_city_id} to {end_city_id}: {path}")
                print(f"Total Distance: {total_distance}")
                print(f'Total Time:  {total_time}')
                print(f'Total Cost: {total_cost + fuel_cost}')   
                print(f'Fuel buy in: {fuel_purchase}')
                print(f'Fuel Cost: {fuel_cost}')

allA_star_cost_output = []
def save_for_cost():
    
    IDs = []
    for id in range(len(A_starcities)):
        IDs.append(A_starcities[id][0])
        
    for i in IDs:
        for j in IDs:
            if i != j:
                start_city_id = int(i)
                end_city_id = int(j)
                path = A_star_path(start_city_id, end_city_id, A_startotalCost, A_stardistances)
                total_distance = calculate_path_distance(path, A_stardistances) 
                total_time = calculate_path_time(path, A_starTime) 
                total_toll = calculate_path_toll(path, A_starTollPrice)
                total_cost = calculate_path_total_cost(path, A_startotalCost)
                fuel_purchase , fuel_cost = minimize_fuel_cost(path, fuel_needed, A_starFuelPrice, F, Sf)

                allA_star_cost_output.append([start_city_id, end_city_id, total_cost + fuel_cost, total_distance, path, total_toll, total_time, fuel_purchase, fuel_cost])
                                
def saveCostSolutionOutput():
    df_SolutionMatrix = pd.DataFrame(allA_star_cost_output)
    # define csv file name
    sol_file_name = 'solution_cost_' + now.strftime("%Y-%m-%d_%H-%M-%S") + '.csv'
    # set column names
    df_SolutionMatrix.columns = ['from', 'to', 'Total Cost', 'pathDistance', 'path', 'toll', 'time', 'City-Amount', 'fuel cost']

    
    if not os.path.exists(sol_directory):
        os.makedirs(sol_directory)
    # writing data frame to a CSV file
    df_SolutionMatrix.to_csv(os.path.join(sol_directory, sol_file_name), sep=';', index=False)
          
def saveSolution():
    readAllAsCSV()
    save_for_cost()
    saveCostSolutionOutput()

def saveDistanceSolution():
    readAllAsCSV()
    save_for_distance()
    saveDistanceSolutionOutput()








