import requests
from bs4 import BeautifulSoup
import pandas as pd
from urllib.parse import urljoin
import time

BASE_URL = "https://www.german-hospital-directory.com"
START_URL = "https://www.german-hospital-directory.com/app/search/state/bavaria"

def get_hospital_data(page_url):
    """Extract hospital data from a single page"""
    try:
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        }
        response = requests.get(page_url, headers=headers)
        response.raise_for_status()
        
        soup = BeautifulSoup(response.text, 'html.parser')
        hospitals = []
        
        # Each hospital is in a div with class 'hospital-item'
        for item in soup.select('div.hospital-item'):
            name = item.select_one('h2.hospital-title a').get_text(strip=True)
            
            # Extract address information
            address_div = item.select_one('div.hospital-address')
            address_lines = [line.get_text(strip=True) for line in address_div.select('p') if line.get_text(strip=True)]
            
            # Some hospitals might have different address formats
            address = ', '.join(address_lines) if address_lines else ''
            
            hospitals.append({
                'name': name,
                'address': address
            })
        
        # Check for next page
        next_page = soup.select_one('li.next a')
        next_url = urljoin(BASE_URL, next_page['href']) if next_page else None
        
        return hospitals, next_url
    
    except Exception as e:
        print(f"Error processing {page_url}: {str(e)}")
        return [], None

def crawl_all_hospitals(start_url):
    """Crawl through all pages of hospital listings"""
    all_hospitals = []
    current_url = start_url
    page_count = 0
    
    while current_url:
        page_count += 1
        print(f"Processing page {page_count}: {current_url}")
        
        hospitals, next_url = get_hospital_data(current_url)
        all_hospitals.extend(hospitals)
        
        current_url = next_url
        
        # Be polite - add delay between requests
        time.sleep(1)
    
    return all_hospitals

def save_to_csv(hospitals, filename='bavaria_hospitals.csv'):
    """Save hospital data to CSV"""
    df = pd.DataFrame(hospitals)
    df.to_csv(filename, index=False, encoding='utf-8-sig')  # utf-8-sig for Excel compatibility
    print(f"Saved {len(df)} hospitals to {filename}")

if __name__ == "__main__":
    print("Starting to crawl hospital data...")
    hospitals = crawl_all_hospitals(START_URL)
    save_to_csv(hospitals)
    print("Crawling completed!")