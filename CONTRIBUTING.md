# ccinput Contribution Guidelines

If you know a quantum chemistry package relatively well and are familiar with Python, you can certainly help the project. We appreciate the efforts that all contributors put into improving `ccinput` and making computational chemistry more accessible and headache-free!

These guidelines will help you contribute and following them will greatly speed up the acceptance of your contribution.

## Features that you can add to the project

Any additional feature is welcome. You can view the [roadmap](https://github.com/orgs/cyllab/projects/1/views/1) to see the planned features; feel free to contribute towards any unassigned task. You can announce that you are starting to work on a feature in the [Discussions tab](https://github.com/cyllab/ccinput/discussions). Feel free to ask for precisions or feedback as well!

General improvements that you can submit:

- Support for new packages (even if rudimentary)
- Support for new features of supported packages
- Correction or improvement of static data (like DFT functionals or basis sets)
- Improvement to the code quality (refactoring, documentation...)
- Addition of tests
- Fix to a bug
- Improvement or extension of the documentation

In addition to this, feature suggestions and bug reports are welcome. Even if you don't have the time to make the changes yourself, reporting the changes to be made will help the project. If you are unuse whether a particular change should be done, feel free to use the [discussions tab](https://github.com/cyllab/ccinput/discussions) to ask for feedback! There are issue templates to guide you through the submission of both feature requests and bug reports. Simply open an issue on Github and choose the relevant template.

## How to contribute
1. Identify the change you want to make (adding a particular feature, resolving an issue, *etc.*).
2. Fork the `ccinput` repository. This will create a copy of the `ccinput` project under your user account (`<username>/ccinput`).
3. Clone your fork on your local machine (`git clone https://github.com/<username>/ccinput.git`).
4. Create a new branch for your modification (`git checkout -b <branch_name>`). Choose a branch name linked with your change.
5. Modify the code to add your feature or fix an issue.
6. Run the automated tests to make sure that everything works as expected (`nosetests` in the main directory or `python setup.py test`). If there are errors, adjust your changes to make the tests pass again.
7. Choose the modified files to be committed (`git add <file>`)
8. Commit the modifications with a short summary of the changes made (`git commit -m "Added the option to [...]"`). Your changes can be split into multiple commits, especially if they are extensive.
9. Push the changes to your fork on Github (`git push origin <branch_name>`)
10. On Github, submit a pull request to the main repository, `cyllab/ccinput`. Add information about the changes made and link to related issues with "#<issue_number>", if any. If you have questions or are unsure about anything, indicate it in the pull request. We can help you fix rough edges in your pull request if you detail what is causing you problems.
11. Wait until we review the pull request and provide feedback or merge the pull request.

Once your pull request is merged, it is officially part of the `ccinput` code and will be included in the next release.

## Code Guidelines

To ensure that the project remains clean and maintainable, it is necessary to establish some guidelines on the code itself. These have been kept to a minimum, but need to be followed. Furthermore, the code style is automatically standardized by [Black](https://github.com/psf/black) when submitting a pull request, and so you do not need to worry about code style.

- Make variable and function names as descriptive yet concise as possible
- Add comments or docstrings to explain what the average contributor will likely not know or understand immediately. However, do not add comments to explain what is obvious to everyone.

**Useful comment:**

	def standardize_xyz(xyz):

	    [...]

	    # Check if the xyz contains the two line header
	    # Remove it if it does
	    try:
		num_atoms = int(arr_xyz[0])
	    except ValueError:
		pass
	    else:
		if num_atoms == len(arr_xyz)-2:
		    arr_xyz = arr_xyz[2:]
		else:
		    raise InvalidXYZ(f"Invalid xyz header: {num_atoms} atoms specified, " +
			    f"but actually contains {len(arr_xyz) - 2} atoms")

	    [...]


**Useless comments:**

	def parse_xyz_from_file(path):
	    """ Given a file path, parses the XYZ structure from the file """

	    if not os.path.isfile(path): # If the file does not exist
		raise InvalidParameter(f"Input file not found: {path}") # Raise an error

	    with open(path) as f: # Open the file
		lines = f.readlines() # Read the XYZ structure from the file

	    return standardize_xyz(lines) # Return the standardized XYZ

Lastly, **please provide at least some tests with your pull request**. If you add a feature, the tests should be enough to verify that the feature works in its most common use cases. If your pull request fixes a bug, provide enough tests to verify that the bug is no longer present after your fix. This ensures that working features keep working and that dead bugs stay dead.

If you are unsure if you have enough tests, you can send a pull request and we will look over it. Any subsequent commits you make to the branch will update the pull request automatically.

Changes to the documentation do not require tests or comments.
