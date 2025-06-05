from flask import Flask, render_template, jsonify, request
from sgp4.api import Satrec, jday
import datetime
import requests
import json
import os
from math import sqrt, atan2, degrees, asin, cos, sin, radians, pi
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import threading
import time
from datetime import datetime, timedelta
import shutil

app = Flask(__name__)

@dataclass
class SatelliteInfo:
    """Satellite metadata and TLE information"""
    norad_id: str
    name: str
    satellite_type: str
    country: str
    launch_date: str
    active: bool
    tle_line1: str
    tle_line2: str
    altitude: float = 0.0
    last_updated: datetime = None

@dataclass
class SpaceCenter:
    """Space center location information"""
    name: str
    country: str
    latitude: float
    longitude: float
    elevation: float = 0.0  # meters above sea level
    min_elevation_angle: float = 10.0  # minimum elevation angle for visibility

class SatelliteTracker:
    """Enhanced satellite tracking system with space center pass predictions"""
    
    def __init__(self):
        self.satellites: Dict[str, SatelliteInfo] = {}
        self.satellite_objects: Dict[str, Satrec] = {}
        self.positions: Dict[str, Dict] = {}
        self.orbital_paths: Dict[str, List[Tuple[float, float]]] = {}
        self.space_centers = self.load_space_centers()
        self.load_satellite_catalog()
        
    def load_space_centers(self) -> Dict[str, SpaceCenter]:
        """Load major space centers around the world"""
        centers = {
            "kennedy": SpaceCenter(
                name="Kennedy Space Center",
                country="USA",
                latitude=28.5721,
                longitude=-80.6480,
                elevation=3.0,
                min_elevation_angle=0.0
            ),
            "baikonur": SpaceCenter(
                name="Baikonur Cosmodrome",
                country="Kazakhstan",
                latitude=45.9200,
                longitude=63.3422,
                elevation=90.0,
                min_elevation_angle=0.0
            ),
            "cape_canaveral": SpaceCenter(
                name="Cape Canaveral Space Force Station",
                country="USA",
                latitude=28.4889,
                longitude=-80.5778,
                elevation=16.0,
                min_elevation_angle=0.0
            ),
            "vandenberg": SpaceCenter(
                name="Vandenberg Space Force Base",
                country="USA",
                latitude=34.7420,
                longitude=-120.5724,
                elevation=110.0,
                min_elevation_angle=0.0
            ),
            "kourou": SpaceCenter(
                name="Guiana Space Centre",
                country="French Guiana",
                latitude=5.2389,
                longitude=-52.7681,
                elevation=50.0,
                min_elevation_angle=0.0
            ),
            "plesetsk": SpaceCenter(
                name="Plesetsk Cosmodrome",
                country="Russia",
                latitude=62.9572,
                longitude=40.5769,
                elevation=118.0,
                min_elevation_angle=0.0
            ),
            "xichang": SpaceCenter(
                name="Xichang Satellite Launch Center",
                country="China",
                latitude=28.2467,
                longitude=102.0264,
                elevation=1825.0,
                min_elevation_angle=0.0
            ),
            "jiuquan": SpaceCenter(
                name="Jiuquan Satellite Launch Center",
                country="China",
                latitude=40.9603,
                longitude=100.2914,
                elevation=1000.0,
                min_elevation_angle=0.0
            ),
            "tanegashima": SpaceCenter(
                name="Tanegashima Space Center",
                country="Japan",
                latitude=30.4008,
                longitude=130.9681,
                elevation=70.0,
                min_elevation_angle=0.0
            ),
            "sriharikota": SpaceCenter(
                name="Satish Dhawan Space Centre",
                country="India",
                latitude=13.7199,
                longitude=80.2304,
                elevation=3.0,
                min_elevation_angle=0.0
            ),
            "palmachim": SpaceCenter(
                name="Palmachim Airbase",
                country="Israel",
                latitude=31.8947,
                longitude=34.6919,
                elevation=45.0,
                min_elevation_angle=0.0
            ),
            "wallops": SpaceCenter(
                name="Wallops Flight Facility",
                country="USA",
                latitude=37.8402,
                longitude=-75.4878,
                elevation=11.0,
                min_elevation_angle=0.0
            ),
            "esa": SpaceCenter(
                name="ESA",
                country="Europe",
                latitude=40.4272,
                longitude=-3.7134,
                elevation=0.0,
                min_elevation_angle=0.0
            ),
            "isro": SpaceCenter(
                name="ISRO",
                country="India",
                latitude=13.7199,
                longitude=80.2305,
                elevation=0.0,
                min_elevation_angle=0.0
            ),
            "nasa": SpaceCenter(
                name="NASA",
                country="USA",
                latitude=28.5721,
                longitude=-80.6480,
                elevation=0.0,
                min_elevation_angle=0.0
            ),
            "svalbard": SpaceCenter(
                name="Svalbard Station",
                country="Norway",
                latitude=78.2306,
                longitude=15.4078,
                elevation=0.0,
                min_elevation_angle=0.0
            )
        }
        
        print(f"📍 Loaded {len(centers)} space centers for tracking")
        return centers
        
    def load_satellite_catalog(self):
        """Load satellite catalog from configuration"""
        catalog = {
            "cartosat2f": {
                "name": "Cartosat-2F",
                "norad_id": "42063",
                "type": "Earth Observation",
                "country": "India",
                "launch_date": "2017-01-12",
                "active": True,
                "tle_file": "cartosat2f.tle"
            },
            "iss": {
                "name": "International Space Station",
                "norad_id": "25544",
                "type": "Space Station",
                "country": "International",
                "launch_date": "1998-11-20",
                "active": True,
                "tle_file": "iss.tle"
            },
            "hubble": {
                "name": "Hubble Space Telescope",
                "norad_id": "20580",
                "type": "Space Telescope",
                "country": "USA",
                "launch_date": "1990-04-24",
                "active": True,
                "tle_file": "hubble.tle"
            },
            "sentinel1a": {
                "name": "Sentinel-1A",
                "norad_id": "39634",
                "type": "Earth Observation",
                "country": "ESA",
                "launch_date": "2014-04-03",
                "active": True,
                "tle_file": "sentinel1a.tle"
            },
            "landsat8": {
                "name": "Landsat 8",
                "norad_id": "39084",
                "type": "Earth Observation",
                "country": "USA",
                "launch_date": "2013-02-11",
                "active": True,
                "tle_file": "landsat8.tle"
            },
            "starlink1": {
                "name": "Starlink-1007",
                "norad_id": "44713",
                "type": "Communication",
                "country": "USA",
                "launch_date": "2019-11-11",
                "active": True,
                "tle_file": "starlink1.tle"
            }
        }
        
        for sat_id, info in catalog.items():
            try:
                tle_path = f"tle/{info['tle_file']}"
                if os.path.exists(tle_path):
                    with open(tle_path, "r") as f:
                        lines = f.readlines()
                    
                    if len(lines) >= 3:
                        satellite_info = SatelliteInfo(
                            norad_id=info['norad_id'],
                            name=info['name'],
                            satellite_type=info['type'],
                            country=info['country'],
                            launch_date=info['launch_date'],
                            active=info['active'],
                            tle_line1=lines[1].strip(),
                            tle_line2=lines[2].strip(),
                            last_updated=datetime.utcnow()
                        )
                        
                        self.satellites[sat_id] = satellite_info
                        self.satellite_objects[sat_id] = Satrec.twoline2rv(
                            satellite_info.tle_line1, 
                            satellite_info.tle_line2
                        )
                        
                        print(f"✅ Loaded satellite: {info['name']}")
                else:
                    print(f"❌ TLE file not found: {tle_path}")
                    # Create sample TLE for demo purposes
                    self.create_sample_tle(sat_id, info)
                    
            except Exception as e:
                print(f"❌ Error loading satellite {sat_id}: {e}")
    
    def create_sample_tle(self, sat_id: str, info: Dict):
        """Create sample TLE data for demonstration"""
        sample_tles = {
            "iss": [
                "1 25544U 98067A   23001.00000000  .00002182  00000-0  40768-4 0  9990",
                "2 25544  51.6461 289.1012 0001354  94.8340 265.3446 15.48919103    10"
            ],
            "hubble": [
                "1 20580U 90037B   23001.00000000  .00000371  00000-0  18543-4 0  9992",
                "2 20580  28.4684 159.4832 0002649  12.8134 347.2785 15.09309617    07"
            ],
            "sentinel1a": [
                "1 39634U 14016A   23001.00000000  .00000023  00000-0  15002-4 0  9993",
                "2 39634  98.1817  61.2976 0001449  85.9467 274.2171 14.59198896    05"
            ],
            "landsat8": [
                "1 39084U 13008A   23001.00000000 -.00000006  00000-0  10189-4 0  9991",
                "2 39084  98.2123  78.4567 0001361  73.8901 286.2701 14.57107527    03"
            ],
            "starlink1": [
                "1 44713U 19074A   23001.00000000  .00001247  00000-0  10034-3 0  9995",
                "2 44713  53.0544 123.4567 0001234  89.1234 270.9876 15.05123456    06"
            ]
        }
        
        if sat_id in sample_tles:
            satellite_info = SatelliteInfo(
                norad_id=info['norad_id'],
                name=info['name'],
                satellite_type=info['type'],
                country=info['country'],
                launch_date=info['launch_date'],
                active=info['active'],
                tle_line1=sample_tles[sat_id][0],
                tle_line2=sample_tles[sat_id][1],
                last_updated=datetime.utcnow()
            )
            
            self.satellites[sat_id] = satellite_info
            self.satellite_objects[sat_id] = Satrec.twoline2rv(
                satellite_info.tle_line1, 
                satellite_info.tle_line2
            )
    
    def calculate_elevation_azimuth(self, sat_lat: float, sat_lon: float, sat_alt: float,
                                  obs_lat: float, obs_lon: float, obs_alt: float = 0) -> Tuple[float, float]:
        """Calculate elevation and azimuth angles from observer to satellite"""
        # Convert to radians
        sat_lat_rad = radians(sat_lat)
        sat_lon_rad = radians(sat_lon)
        obs_lat_rad = radians(obs_lat)
        obs_lon_rad = radians(obs_lon)
        
        # Earth radius in km
        R = 6371.0
        
        # Convert satellite position to Cartesian coordinates
        sat_x = (R + sat_alt) * cos(sat_lat_rad) * cos(sat_lon_rad)
        sat_y = (R + sat_alt) * cos(sat_lat_rad) * sin(sat_lon_rad)
        sat_z = (R + sat_alt) * sin(sat_lat_rad)
        
        # Convert observer position to Cartesian coordinates
        obs_x = (R + obs_alt / 1000) * cos(obs_lat_rad) * cos(obs_lon_rad)
        obs_y = (R + obs_alt / 1000) * cos(obs_lat_rad) * sin(obs_lon_rad)
        obs_z = (R + obs_alt / 1000) * sin(obs_lat_rad)
        
        # Vector from observer to satellite
        dx = sat_x - obs_x
        dy = sat_y - obs_y
        dz = sat_z - obs_z
        
        # Distance to satellite
        distance = sqrt(dx**2 + dy**2 + dz**2)
        
        # Calculate elevation angle
        dot_product = dx * obs_x + dy * obs_y + dz * obs_z
        elevation = degrees(asin(dot_product / (distance * (R + obs_alt / 1000))))
        
        # Calculate azimuth angle (simplified)
        azimuth = degrees(atan2(dy, dx))
        if azimuth < 0:
            azimuth += 360
        
        return elevation, azimuth
    
    def predict_next_passes(self, sat_id: str, center_id: str, 
                           duration_hours: int = 87600) -> List[Dict]:
        """Predict next passes of a satellite over a space center"""
        if sat_id not in self.satellite_objects or center_id not in self.space_centers:
            return []
        
        satellite = self.satellite_objects[sat_id]
        center = self.space_centers[center_id]
        passes = []
        
        now = datetime.utcnow()
        current_pass = None
        
        # Check every minute for the specified duration
        for i in range(duration_hours * 60):
            check_time = now + timedelta(minutes=i)
            jd, fr = jday(check_time.year, check_time.month, check_time.day,
                         check_time.hour, check_time.minute, check_time.second)
            
            e, r, v = satellite.sgp4(jd, fr)
            
            if e == 0:  # No error
                x, y, z = r
                R = sqrt(x**2 + y**2 + z**2)
                sat_lat = degrees(asin(z / R))
                sat_lon = degrees(atan2(y, x))
                sat_alt = R - 6371.0
                
                # Calculate elevation and azimuth
                elevation, azimuth = self.calculate_elevation_azimuth(
                    sat_lat, sat_lon, sat_alt,
                    center.latitude, center.longitude, center.elevation
                )
                
                # Check if satellite is visible (above minimum elevation)
                if elevation >= center.min_elevation_angle:
                    if current_pass is None:
                        # Start of a new pass
                        current_pass = {
                            'start_time': check_time,
                            'max_elevation': elevation,
                            'max_elevation_time': check_time,
                            'max_elevation_azimuth': azimuth,
                            'duration_minutes': 1
                        }
                    else:
                        # Continue current pass
                        current_pass['duration_minutes'] += 1
                        if elevation > current_pass['max_elevation']:
                            current_pass['max_elevation'] = elevation
                            current_pass['max_elevation_time'] = check_time
                            current_pass['max_elevation_azimuth'] = azimuth
                else:
                    if current_pass is not None:
                        # End of current pass
                        current_pass['end_time'] = check_time - timedelta(minutes=1)
                        passes.append({
                            'satellite_id': sat_id,
                            'satellite_name': self.satellites[sat_id].name,
                            'space_center': center.name,
                            'start_time': current_pass['start_time'].strftime("%Y-%m-%d %H:%M:%S UTC"),
                            'end_time': current_pass['end_time'].strftime("%Y-%m-%d %H:%M:%S UTC"),
                            'max_elevation_time': current_pass['max_elevation_time'].strftime("%Y-%m-%d %H:%M:%S UTC"),
                            'max_elevation_degrees': round(current_pass['max_elevation'], 2),
                            'max_elevation_azimuth': round(current_pass['max_elevation_azimuth'], 2),
                            'duration_minutes': current_pass['duration_minutes'],
                            'visibility': 'Excellent' if current_pass['max_elevation'] > 60 else 
                                        'Good' if current_pass['max_elevation'] > 30 else 'Fair'
                        })
                        current_pass = None
        
        # Handle case where pass extends beyond our check period
        if current_pass is not None:
            current_pass['end_time'] = now + timedelta(hours=duration_hours)
            passes.append({
                'satellite_id': sat_id,
                'satellite_name': self.satellites[sat_id].name,
                'space_center': center.name,
                'start_time': current_pass['start_time'].strftime("%Y-%m-%d %H:%M:%S UTC"),
                'end_time': current_pass['end_time'].strftime("%Y-%m-%d %H:%M:%S UTC"),
                'max_elevation_time': current_pass['max_elevation_time'].strftime("%Y-%m-%d %H:%M:%S UTC"),
                'max_elevation_degrees': round(current_pass['max_elevation'], 2),
                'max_elevation_azimuth': round(current_pass['max_elevation_azimuth'], 2),
                'duration_minutes': current_pass['duration_minutes'],
                'visibility': 'Excellent' if current_pass['max_elevation'] > 60 else 
                            'Good' if current_pass['max_elevation'] > 30 else 'Fair'
            })
        
        return passes
    
    def get_all_upcoming_passes(self, duration_hours: int = 87600) -> Dict[str, List[Dict]]:
        """Get all upcoming passes for all satellites over all space centers"""
        all_passes = {}
        
        for sat_id in self.satellites.keys():
            sat_passes = []
            for center_id in self.space_centers.keys():
                passes = self.predict_next_passes(sat_id, center_id, duration_hours)
                sat_passes.extend(passes)
            
            # Sort passes by start time
            sat_passes.sort(key=lambda x: x['start_time'])
            all_passes[sat_id] = sat_passes
        
        return all_passes
    
    def get_next_pass_summary(self) -> List[Dict]:
        """Get summary of next passes for all satellites"""
        summary = []
        now = datetime.utcnow()
        
        for sat_id, sat_info in self.satellites.items():
            next_passes = []
            
            for center_id, center in self.space_centers.items():
                passes = self.predict_next_passes(sat_id, center_id, 87600)
                if passes:
                    next_passes.append(passes[0])  # Get the first (earliest) pass
            
            if next_passes:
                # Sort by start time and get the earliest
                next_passes.sort(key=lambda x: x['start_time'])
                earliest_pass = next_passes[0]
                
                summary.append({
                    'satellite_id': sat_id,
                    'satellite_name': sat_info.name,
                    'satellite_type': sat_info.satellite_type,
                    'next_pass': earliest_pass,
                    'time_until_pass': self.calculate_time_until(earliest_pass['start_time'])
                })
        
        return summary
    
    def calculate_time_until(self, future_time_str: str) -> str:
        """Calculate human-readable time until a future event"""
        future_time = datetime.strptime(future_time_str, "%Y-%m-%d %H:%M:%S UTC")
        now = datetime.utcnow()
        delta = future_time - now
        
        if delta.total_seconds() < 0:
            return "Past"
        
        hours = int(delta.total_seconds() // 3600)
        minutes = int((delta.total_seconds() % 3600) // 60)
        
        if hours > 0:
            return f"{hours}h {minutes}m"
        else:
            return f"{minutes}m"
    
    def get_satellite_position(self, sat_id: str) -> Optional[Dict]:
        """Get current position of a specific satellite"""
        if sat_id not in self.satellite_objects:
            return None
        
        satellite = self.satellite_objects[sat_id]
        now = datetime.utcnow()
        jd, fr = jday(now.year, now.month, now.day, now.hour, now.minute, now.second)
        
        e, r, v = satellite.sgp4(jd, fr)
        
        if e != 0:
            print(f"❌ Error computing position for {sat_id}: {e}")
            return None
        
        x, y, z = r
        vx, vy, vz = v
        
        # Calculate latitude, longitude, altitude
        R = sqrt(x**2 + y**2 + z**2)
        lat = degrees(asin(z / R))
        lon = degrees(atan2(y, x))
        alt = R - 6371.0  # Earth radius in km
        
        # Calculate velocity
        velocity = sqrt(vx**2 + vy**2 + vz**2)
        
        satellite_info = self.satellites[sat_id]
        
        position_data = {
            'id': sat_id,
            'name': satellite_info.name,
            'type': satellite_info.satellite_type,
            'country': satellite_info.country,
            'active': satellite_info.active,
            'latitude': round(lat, 6),
            'longitude': round(lon, 6),
            'altitude': round(alt, 2),
            'velocity': round(velocity, 2),
            'timestamp': now.strftime("%Y-%m-%d %H:%M:%S UTC"),
            'launch_date': satellite_info.launch_date,
            'norad_id': satellite_info.norad_id
        }
        
        self.positions[sat_id] = position_data
        return position_data
    
    def get_all_positions(self) -> Dict[str, Dict]:
        """Get positions of all satellites as {sat_id: {lat, lon, altitude}} using ECI conversion"""
        results = {}
        now = datetime.utcnow()
        for sat_id in self.satellites.keys():
            if sat_id in self.satellite_objects:
                satellite = self.satellite_objects[sat_id]
                jd, fr = jday(now.year, now.month, now.day, now.hour, now.minute, now.second)
                e, r, v = satellite.sgp4(jd, fr)
                if e == 0:
                    x, y, z = r
                    R = sqrt(x**2 + y**2 + z**2)
                    lat = degrees(asin(z / R))
                    lon = degrees(atan2(y, x))
                    alt = R - 6371.0  # Earth radius in km
                    results[sat_id] = {"lat": lat, "lon": lon, "altitude": round(alt, 2)}
        return results
    
    def generate_orbital_path(self, sat_id: str, duration_hours: int = 2, 
                           points: int = 50) -> List[Tuple[float, float]]:
        """Generate predicted orbital path for a satellite"""
        if sat_id not in self.satellite_objects:
            return []
        
        satellite = self.satellite_objects[sat_id]
        now = datetime.utcnow()
        path_points = []
        
        for i in range(points):
            # Calculate time offset
            time_offset = timedelta(hours=(duration_hours * i / points))
            future_time = now + time_offset
            
            jd, fr = jday(future_time.year, future_time.month, future_time.day,
                         future_time.hour, future_time.minute, future_time.second)
            
            e, r, v = satellite.sgp4(jd, fr)
            
            if e == 0:  # No error
                x, y, z = r
                R = sqrt(x**2 + y**2 + z**2)
                lat = degrees(asin(z / R))
                lon = degrees(atan2(y, x))
                path_points.append((lat, lon))
        
        self.orbital_paths[sat_id] = path_points
        return path_points
    
    def calculate_distance(self, lat1: float, lon1: float, 
                          lat2: float, lon2: float) -> float:
        """Calculate distance between two points on Earth"""
        R = 6371.0  # Earth radius in km
        
        lat1_rad = radians(lat1)
        lon1_rad = radians(lon1)
        lat2_rad = radians(lat2)
        lon2_rad = radians(lon2)
        
        dlat = lat2_rad - lat1_rad
        dlon = lon2_rad - lon1_rad
        
        a = sin(dlat/2)**2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)**2
        c = 2 * atan2(sqrt(a), sqrt(1-a))
        
        return R * c
    
    def update_tle_data(self):
        """Update TLE data from online sources"""
        print("🔄 Updating TLE data from CelesTrak...")
        
        # URLs for different satellite categories
        tle_sources = {
            "stations": "https://celestrak.com/NORAD/elements/stations.txt",
            "visual": "https://celestrak.com/NORAD/elements/visual.txt",
            "science": "https://celestrak.com/NORAD/elements/science.txt",
            "earth-resources": "https://celestrak.com/NORAD/elements/resource.txt"
        }
        
        try:
            for category, url in tle_sources.items():
                response = requests.get(url, timeout=30)
                if response.status_code == 200:
                    # Parse TLE data and update relevant satellites
                    tle_lines = response.text.strip().split('\n')
                    
                    for i in range(0, len(tle_lines), 3):
                        if i + 2 < len(tle_lines):
                            name = tle_lines[i].strip()
                            line1 = tle_lines[i + 1].strip()
                            line2 = tle_lines[i + 2].strip()
                            
                            # Update if we're tracking this satellite
                            for sat_id, sat_info in self.satellites.items():
                                if sat_info.name.upper() in name.upper():
                                    sat_info.tle_line1 = line1
                                    sat_info.tle_line2 = line2
                                    sat_info.last_updated = datetime.utcnow()
                                    
                                    # Update SGP4 object
                                    self.satellite_objects[sat_id] = Satrec.twoline2rv(line1, line2)
                                    print(f"✅ Updated TLE for {sat_info.name}")
            
            print("✅ TLE update completed")
            
        except Exception as e:
            print(f"❌ Error updating TLE data: {e}")

    def eci_to_latlon(self, r: Tuple[float, float, float], dt: datetime):
        x, y, z = r
        # Convert to latitude, longitude, altitude (basic estimate)
        r_mag = sqrt(x**2 + y**2 + z**2)
        lat = asin(z / r_mag)
        lon = atan2(y, x)
        return degrees(lat), degrees(lon)

def download_global_tle(tle_path='tle/all_active.txt', url='https://celestrak.org/NORAD/elements/gp.php?GROUP=active&FORMAT=tle'):
    """Download the global TLE file if not present or older than 1 day."""
    try:
        if not os.path.exists('tle'):
            os.makedirs('tle')
        if os.path.exists(tle_path):
            mtime = os.path.getmtime(tle_path)
            age_hours = (time.time() - mtime) / 3600
            if age_hours < 24:
                return  # Already up to date
        print(f"🌐 Downloading global TLE from {url} ...")
        r = requests.get(url, timeout=60)
        if r.status_code == 200:
            with open(tle_path, 'w', encoding='utf-8') as f:
                f.write(r.text)
            print(f"✅ Downloaded and saved {tle_path}")
        else:
            print(f"❌ Failed to download TLE: HTTP {r.status_code}")
    except Exception as e:
        print(f"❌ Error downloading global TLE: {e}")

@app.route('/')
def home():
    cartosat_position = tracker.get_satellite_position('cartosat2f')
    position = None
    if cartosat_position:
        position = [cartosat_position['latitude'], cartosat_position['longitude']]
    return render_template("index.html", position=position)

@app.route('/api/satellites')
def api_satellites():
    """Get list of all tracked satellites"""
    satellites = []
    for sat_id, sat_info in tracker.satellites.items():
        satellites.append({
            'id': sat_id,
            'name': sat_info.name,
            'type': sat_info.satellite_type,
            'country': sat_info.country,
            'active': sat_info.active,
            'launch_date': sat_info.launch_date,
            'norad_id': sat_info.norad_id,
            'last_updated': sat_info.last_updated.strftime("%Y-%m-%d %H:%M:%S UTC") if sat_info.last_updated else None
        })
    
    return jsonify({'satellites': satellites, 'count': len(satellites)})

@app.route('/api/space-centers')
def api_space_centers():
    """Get list of all space centers"""
    centers = []
    for center_id, center in tracker.space_centers.items():
        centers.append({
            'id': center_id,
            'name': center.name,
            'country': center.country,
            'latitude': center.latitude,
            'longitude': center.longitude,
            'elevation': center.elevation,
            'min_elevation_angle': center.min_elevation_angle
        })
    
    return jsonify({'space_centers': centers, 'count': len(centers)})

@app.route('/api/positions')
def api_positions():
    """Get current positions of all satellites"""
    positions = tracker.get_all_positions()
    
    return jsonify({
        'positions': positions,
        'count': len(positions),
        'timestamp': datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    })

@app.route('/api/position/<sat_id>')
def api_position(sat_id):
    """Get position of a specific satellite"""
    position = tracker.get_satellite_position(sat_id)
    
    if position:
        return jsonify(position)
    else:
        return jsonify({'error': f'Satellite {sat_id} not found or error in calculation'}), 404

@app.route('/api/next-passes')
def api_next_passes():
    """Get summary of next passes for all satellites"""
    summary = tracker.get_next_pass_summary()
    
    return jsonify({
        'next_passes': summary,
        'count': len(summary),
        'generated_at': datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    })

@app.route('/api/passes/<sat_id>/<center_id>')
def api_specific_passes(sat_id, center_id):
    """Get passes of a specific satellite over a specific space center"""
    duration = request.args.get('duration', default=24, type=int)
    passes = tracker.predict_next_passes(sat_id, center_id, duration)
    return jsonify({
        'satellite_id': sat_id,
        'center_id': center_id,
        'passes': passes,
        'count': len(passes),
        'generated_at': datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S UTC")
    })

@app.route('/api/update-tle', methods=['POST', 'GET'])
def api_update_tle():
    tracker.update_tle_data()
    return jsonify({'status': 'TLE update triggered'})

@app.route('/api/control_station_passes')
def control_station_passes():
    station_param = request.args.get('station', '').lower()
    now = datetime.utcnow()
    station_results = {}

    for station in CONTROL_STATIONS:
        station_name = station['name']
        if station_param and station_param != station_name.lower():
            continue  # 🔥 Only process the requested station if a filter is given

        center_obj = None
        for center in tracker.space_centers.values():
            if center.name.lower() == station['name'].lower():
                center_obj = center
                break

        passes = []
        for sat_id, sat_info in tracker.satellites.items():
            center_id = None
            for cid, c in tracker.space_centers.items():
                if c.name.lower().startswith(station_name.lower()):
                    center_id = cid
                    break
            if center_id:
                sat_passes = tracker.predict_next_passes(sat_id, center_id, 24)
                if sat_passes:
                    p = sat_passes[0]
                    passes.append({
                        'satellite': sat_info.name,
                        'satellite_id': sat_id,
                        'start_time': p['start_time'],
                        'max_elevation_time': p['max_elevation_time'],
                        'max_elevation_degrees': p['max_elevation_degrees'],
                        'visibility': p['visibility'],
                    })

        passes.sort(key=lambda x: x['start_time'])
        next_passes = passes[:3]  # 🎯 Show up to 3 upcoming passes
        past_passes = []  # You can enhance this later

        station_results[station_name.lower()] = {
            'next_passes': next_passes,
            'past_passes': past_passes
        }
        print(f"\n📡 Station Requested: {station_param}")
        print(json.dumps(station_results, indent=2))

    return jsonify({'stations': station_results})

@app.route('/api/search-satellites')
def api_search_satellites():
    """Search all satellites in the global TLE file by name or NORAD ID."""
    query = request.args.get('q', '').strip().lower()
    if not query:
        return jsonify({'error': 'No query provided'}), 400
    tle_path = 'tle/all_active.txt'
    download_global_tle(tle_path)
    results = []
    try:
        with open(tle_path, 'r', encoding='utf-8') as f:
            # Read and filter out blank lines
            lines = [line.strip() for line in f if line.strip()]
        i = 0
        while i < len(lines) - 2:
            name = lines[i]
            line1 = lines[i+1]
            line2 = lines[i+2]
            # Ensure correct TLE format
            if line1.startswith('1 ') and line2.startswith('2 '):
                norad_id = line1.split()[1] if len(line1.split()) > 1 else ''
                norad_id_stripped = norad_id.rstrip('U').lstrip('0')
                query_stripped = query.rstrip('u').lstrip('0')
                # Match by name (case-insensitive, stripped)
                name_match = query in name.lower()
                # Match by NORAD ID (with or without trailing U, leading zeros)
                norad_match = (
                    query == norad_id.lower() or
                    query == norad_id_stripped.lower() or
                    query_stripped == norad_id.lower() or
                    query_stripped == norad_id_stripped.lower()
                )
                if name_match or norad_match:
                    results.append({
                        'name': name,
                        'norad_id': norad_id,
                        'tle_line1': line1,
                        'tle_line2': line2
                    })
                    if len(results) >= 20:
                        break  # Limit results for performance
                i += 3
            else:
                i += 1  # Skip to next line if not a valid TLE group
        return jsonify({'results': results, 'count': len(results)})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

@app.route('/api/position-by-tle', methods=['POST'])
def api_position_by_tle():
    """Compute current position from TLE lines (for searched satellites)."""
    data = request.get_json()
    if not data or 'line1' not in data or 'line2' not in data:
        return jsonify({'error': 'Missing TLE lines'}), 400
    line1 = data['line1']
    line2 = data['line2']
    try:
        sat = Satrec.twoline2rv(line1, line2)
        now = datetime.utcnow()
        jd, fr = jday(now.year, now.month, now.day, now.hour, now.minute, now.second)
        e, r, v = sat.sgp4(jd, fr)
        if e != 0:
            return jsonify({'error': f'SGP4 error code: {e}'}), 400
        x, y, z = r
        R = sqrt(x**2 + y**2 + z**2)
        lat = degrees(asin(z / R))
        lon = degrees(atan2(y, x))
        alt = R - 6371.0  # Earth radius in km
        return jsonify({
            'lat': lat,
            'lon': lon,
            'altitude': alt
        })
    except Exception as ex:
        return jsonify({'error': str(ex)}), 500

# Initialize the tracker
tracker = SatelliteTracker()

# Background TLE update (runs every 6 hours)
def background_tle_update():
    while True:
        time.sleep(6 * 3600)  # 6 hours
        tracker.update_tle_data()

tle_update_thread = threading.Thread(target=background_tle_update, daemon=True)
tle_update_thread.start()

# Example control stations (add more as needed)
CONTROL_STATIONS = [
    {
        'name': 'ISRO',
        'lat': 13.0649,
        'lon': 77.6336
    },
    {
        'name': 'NASA',
        'lat': 28.5721,
        'lon': -80.6480
    },
    {
        'name': 'ESA',
        'lat': 40.4272,
        'lon': -3.7134
    },
    {
        'name': 'JAXA',
        'lat': 30.4008,
        'lon': 130.9681
    },
    {
        'name': 'Roscosmos',
        'lat': 62.9572,
        'lon': 40.5769
    },
    {
        'name': 'CNSA',
        'lat': 28.2467,
        'lon': 102.0264
    },
    {
        'name': 'ISRAEL',
        'lat': 31.8947,
        'lon': 34.6919
    },
    {
        'name': 'Arianespace',
        'lat': 5.2389,
        'lon': -52.7681
    },
    {
        'name': 'Wallops',
        'lat': 37.8402,
        'lon': -75.4878
    }
]

def predict_passes(station, satellites, now=None, horizon_deg=10, max_results=3):
    # Dummy implementation: finds satellites currently above the horizon and next 3 closest by distance
    # In a real app, use SGP4 or pyorbital to predict actual pass times
    if now is None:
        now = datetime.utcnow()
    results = []
    past = []
    for sat in satellites:
        # Calculate distance from station to satellite subpoint
        dlat = sat['lat'] - station['lat']
        dlon = sat['lon'] - station['lon']
        dist = (dlat**2 + dlon**2)**0.5
        # Dummy: if within 10 deg, it's overhead now
        if abs(dlat) < 10 and abs(dlon) < 10:
            results.append({
                'satellite': sat['name'],
                'time': now.strftime('%H:%M UTC'),
                'status': 'now'
            })
        else:
            # Predict next pass in X minutes (dummy: use distance as minutes)
            pass_time = now + timedelta(minutes=int(dist*6))
            results.append({
                'satellite': sat['name'],
                'time': pass_time.strftime('%H:%M UTC'),
                'status': 'upcoming'
            })
    # Sort by soonest time, only keep 3
    results.sort(key=lambda x: x['time'])
    upcoming = [r for r in results if r['status'] == 'upcoming'][:max_results]
    now_sats = [r for r in results if r['status'] == 'now']
    # Dummy: move any satellites with time in the past to 'past'
    for r in results:
        t = datetime.strptime(r['time'], '%H:%M UTC')
        if t < now:
            past.append(r)
    return {
        'now': now_sats,
        'upcoming': upcoming,
        'past': past
    }

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)