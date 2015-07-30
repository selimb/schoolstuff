def maketable(input_latex_string,filename):
		"""Modifies input latex string and appends it to filename
		This is the general usecase"""

		inside = input_latex_string.split('\endlastfoot')[1]
		values = inside+'\end{table}'
		lines = values.splitlines()
		for i,line in enumerate(lines):
				if '- T' in line or 'Top' in line: lines[i]+='\\midrule'
		values = '\n'.join(lines)
		pre = open(filename,'r').read().split('\endlastfoot')[0]
		final = pre + '\\endlastfoot' + values
		open(filename,'w').write(final)

def make_stress_table(input_latex_string,filename, is_reversed = False):
	"""Modifies input latex STRESS string and appends it to filename.
	Adds phantom columns where necessary."""
	inside = input_latex_string.split('\endlastfoot')[1]
	values = inside+'\end{table}'
	lines = values.splitlines()
	for i,line in enumerate(lines):
			if '&' not in line:
				continue
			d = '&'
			amber_sep = line.split('&')
			amber_sep[0] +=  '&'
			amber_sep[0:-1] = [r + '&' for r in amber_sep[0:-1]]
			lines[i] = ''.join(amber_sep)
			if not is_reversed:
				if '- T' in line or 'Top' in line: lines[i]+='\\midrule'
			elif is_reversed:
				if '- B' in line or 'Bot' in line: lines[i]+='\\midrule'


	values = '\n'.join(lines)
	pre = open(filename,'r').read().split('\endlastfoot')[0]
	final = pre + '\\endlastfoot' + values
	open(filename,'w').write(final)

def make_safety_table(input_latex_string,filename,is_reversed = False):
	"""Modifies input latex SAFETY string and appends it to filename.
	Adds phantom columns where necessary."""
	inside = input_latex_string.split('\endlastfoot')[1]
	values = inside+'\end{table}'
	lines = values.splitlines()
	for i,line in enumerate(lines):
			if '&' not in line:
				continue
			d = '&'
			amber_sep = line.split('&')
			for j in [0,5,7]:
				amber_sep[j] += '&'
			amber_sep[0:-1] = [r + '&' for r in amber_sep[0:-1]]
			lines[i] = ''.join(amber_sep)
			if not is_reversed:
				if '- T' in line or 'Top' in line: lines[i]+='\\midrule'
			elif is_reversed:
				if '- B' in line or 'Bot' in line: lines[i]+='\\midrule'

	values = '\n'.join(lines)
	pre = open(filename,'r').read().split('\endlastfoot')[0]
	final = pre + '\\endlastfoot' + values
	open(filename,'w').write(final)


if __name__ == '__main__':
		maketable(open('docs/ass4/thetable.tex').read(),'docs/ass4/thetable.tex')
