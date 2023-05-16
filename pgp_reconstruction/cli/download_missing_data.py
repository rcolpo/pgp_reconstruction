#pgp_reconstruction needs very big files, which can not be included in the pip installed version. After installation, on the first time the program in ran, this file downloads the files from a google drive folder.
from pgp_reconstruction import project_dir
import os
import requests
import platform
from bs4 import BeautifulSoup
import zipfile
from datetime import datetime
import shutil
import time
import re

def parse_file_size(file_size):
	size = float(re.findall(r"\d+\.\d+", file_size)[0])
	if "KB" in file_size:
		size *= 1024
	elif "MB" in file_size:
		size *= 1024 ** 2
	elif "GB" in file_size:
		size *= 1024 ** 3
	return size
	
def get_file_size(local_path):
	return os.path.getsize(local_path) if os.path.exists(local_path) else 0
	
def download_file(file_url, local_path, expected_size):

	attempts = 3
	while attempts:
		downloaded_size = 0
		lastPercert = -1
		lastPercertTime = 0
		try:
			response = requests.get(file_url, stream=True)
			with open(local_path, 'wb') as file:
				for chunk in response.iter_content(chunk_size=1024): 
					if chunk:
						file.write(chunk)
						downloaded_size += len(chunk)
						percent = round((downloaded_size / expected_size) * 100, 2)
						if int(percent) <= lastPercert: continue
						if not lastPercertTime: pass
						else:
							delta = datetime.now() - lastPercertTime
							if delta.total_seconds() < 2: continue
						lastPercert = int(percent)
						lastPercertTime = datetime.now()
						print(f"Downloaded: {percent}%")
		except: pass
		
		if downloaded_size < expected_size:
			print("Incomplete download, retrying...")
			time.sleep(5)
			attempts -= 1
		else:
			print("Download complete.")
			break

	return downloaded_size


def download_missing_files():
	#check if newer files exist in the sftp server, comparing with local copy. download files from sftp server if local copy is not found or if remote copy is newer.  
	
	base_url = "https://files.ufz.de/~umb-pgp_reconstruction-01/index.html"
		
	response = requests.get(base_url)
	html_doc = response.text

	soup = BeautifulSoup(html_doc, 'html.parser')
	
	#filePath = 'C:/Users/colpoama/AppData/Roaming/Python/Python37/site-packages/pgp_reconstruction/cli'
	filePath = os.path.realpath(__file__)
	local_base_dir = os.path.dirname(filePath).split('site-packages')[0] + 'site-packages' # go to where packages are saved
	while local_base_dir[-1] == '/' or local_base_dir[-1] == '\\': local_base_dir = local_base_dir[:-1]
	
	
	minPathMissing = 0
	for h2 in soup.find_all('h2'):
		local_dir_ToSaveFiles = os.path.join(local_base_dir + h2.text.replace("Files to be saved on: [...]PythonX\Lib\site-packages", ''))
		table = h2.find_next_sibling('table')
		for row in table.find_all('tr')[1:]:  # skip the header row
		
			cols = row.find_all('td')
			file_url = cols[0].a['href']
			file_name = os.path.basename(file_url)
			local_file = os.path.join(local_dir_ToSaveFiles, file_name)
			last_modified = datetime.strptime(cols[1].text, "%a %b %d %H:%M:%S %Y")
			file_size = cols[2].text # in KB, MB, or GB
			
			expected_size = parse_file_size(file_size)
			
			
			local_file_size = get_file_size(local_file)
			
			if file_name == 'MinPath_master.zip':
				minPath_folder = os.path.join(local_dir_ToSaveFiles, 'MinPath_master')
				
				if not os.path.exists(minPath_folder) or  datetime.fromtimestamp(os.path.getmtime(minPath_folder)) < last_modified:
				
					if os.path.exists(minPath_folder): shutil.rmtree(minPath_folder)
					else: minPathMissing = 1
					
					print(f"Downloading file: {file_url}")
					download_file(file_url, local_file, expected_size)
					
					#unzip minPath and remove zip file
					with zipfile.ZipFile(os.path.join(project_dir, 'dependencies', 'MinPath_master.zip'), 'r') as zip_ref:
						zip_ref.extractall(os.path.join(project_dir, 'dependencies'))
						
					try: os.remove(local_file)
					except: pass
					
			elif file_name == 'cplex_solver.py':
				print(f"Downloading file: {file_url}")
				download_file(file_url, local_file, expected_size)
						
			else:
				if local_file_size < expected_size or datetime.fromtimestamp(os.path.getmtime(local_file)) < last_modified:
				
					print(f"Downloading file: {file_url}")
					downloaded_size = download_file(file_url, local_file, expected_size)


	#check if prodigal is available. If not, download official release
	
	dependencies_folder = os.path.join(local_base_dir, 'pgp_reconstruction', 'dependencies')
	dependencies_files = os.listdir(dependencies_folder)
	
	prodigal = False
	for file in dependencies_files:
		if 'prodigal' in file: prodigal = True
	if prodigal == False:	
		if platform.system() == 'Windows':
			prodigalURL = 'https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.windows.exe'
		if platform.system() == 'Darwin':
			prodigalURL = 'https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.osx.10.9.5'
		if platform.system() == 'Linux':
			prodigalURL = 'https://github.com/hyattpd/Prodigal/releases/download/v2.6.3/prodigal.linux'

		response = requests.get(prodigalURL, stream=True)
		response.raise_for_status()
		save_path = os.path.join(project_dir, 'dependencies', prodigalURL.split('/')[-1])
		with open(save_path, 'wb') as file:
			for chunk in response.iter_content(chunk_size=8192):
				file.write(chunk)
	
	print('\nSync complete.\n')
	
	return minPathMissing
