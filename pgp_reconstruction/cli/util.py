import os

def saveProgressFile(progressValue, outputfolder):
	#Crete a file with the task's approximate progress. It is useful only for the web application. So the web application can use the same package version as in GitHub.
	
	progress_path = os.path.join(outputfolder, 'progress.txt')
	
	with open(progress_path, "w") as f:
		f.write(str(progressValue))
	
	return