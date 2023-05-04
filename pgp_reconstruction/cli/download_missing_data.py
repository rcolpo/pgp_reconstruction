#pgp_reconstruction needs very big files, which can not be included in the pip installed version. After installation, on the first time the program in ran, this file downloads the files from a google drive folder.

from pgp_reconstruction import project_dir
import os
import requests
import platform
from bs4 import BeautifulSoup
import zipfile

def get_google_drive_files(folder_url):
	response = requests.get(folder_url)
	if response.status_code == 200:
		soup = BeautifulSoup(response.text, 'lxml')
		file_elements = soup.find_all('div', class_='WYuW0e Ss7qXc')

		files = []
		for element in file_elements:
			file_name_element = element.find('div', class_='KL4NAf')
			if file_name_element:
				file_name = file_name_element['data-tooltip'].split(': ')[-1]
			else:
				continue
			file_id = element['data-id']
			files.append((file_name, file_id))
		return files
	else:
		print(f'Error: Unable to fetch files from Google Drive. Status code: {response.status_code}')
		return []


def download_file(file_name, file_id, download_folder):
	url = f'https://drive.google.com/uc?export=download&id={file_id}'
	response = requests.get(url)

	if response.status_code == 200:
		with open(os.path.join(download_folder, file_name), 'wb') as f:
			f.write(response.content)
	else:
		print(f'Error: Unable to download file {file_name}. Status code: {response.status_code}')


def download_missing_files():
	#download databases
	data_folder_url = 'https://drive.google.com/drive/u/1/folders/1u1m6mCH8s2gXrWUSxEmk4vMNo4beaRxX'
	google_drive_files = get_google_drive_files(data_folder_url)
	
	data_folder = os.path.join(project_dir, 'data/generated')
	data_files = os.listdir(data_folder)

	for file_name, file_id in google_drive_files:
		if file_name not in data_files:
			print(f'Downloading {file_name}...')
			download_file(file_name, file_id, data_folder)
			
	#download dependencies
	dependencies_folder_url = 'https://drive.google.com/drive/u/1/folders/1oBhyr06whJl_Np0AJ8yPcoeaYrXbeFj6'
	google_drive_files = get_google_drive_files(dependencies_folder_url)
	
	dependencies_folder = os.path.join(project_dir, 'dependencies')
	dependencies_files = os.listdir(dependencies_folder)

	minPathMissing = 1
	if 'MinPath_master' not in dependencies_files:
		for file_name, file_id in google_drive_files:
			continue
			if file_name == 'MinPath_master.zip':
				#Download modified version of minPath
				download_file(file_name, file_id, dependencies_folder)
				
				#unzip minPath and remove zip file
				with zipfile.ZipFile(os.path.join(project_dir, 'dependencies', 'MinPath_master.zip'), 'r') as zip_ref:
					zip_ref.extractall(os.path.join(project_dir, 'dependencies'))
					
				os.remove(os.path.join(project_dir, 'dependencies', 'MinPath_master.zip'))
	else: minPathMissing = 0
		
	
	#check if prodigal is available. If not, download official release
	prodigal = False
	for file in dependencies_folder:
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
	
	print('Sync complete.')
	
	return minPathMissing
	


