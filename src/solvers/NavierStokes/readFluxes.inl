template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::readFluxes()
{
	PetscErrorCode    ierr;
	PetscViewer       viewer;
	std::stringstream ss;
	std::string       savePointDir, fileName;

	ierr = PetscPrintf(PETSC_COMM_WORLD, "Restarting from time step %d.\n", timeStep); CHKERRQ(ierr);

	ss << caseFolder << "/" << std::setfill('0') << std::setw(7) << timeStep;
	savePointDir = ss.str();

	ss.str("");
	ss.clear();
	ss << savePointDir << "/q.dat";
	fileName = ss.str();
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, fileName.c_str(), FILE_MODE_READ, &viewer); CHKERRQ(ierr);
	ierr = VecLoad(q, viewer); CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&viewer); CHKERRQ(ierr);
	
	return 0;
}