template <PetscInt dim>
PetscErrorCode NavierStokesSolver<dim>::createDMs()
{
	PetscErrorCode    ierr;
	DMDABoundaryType  bx, by, bz;
	
	// set boundary types
	bx = (flowDesc->bc[0][XPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	by = (flowDesc->bc[0][YPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;
	if(dim==3) bz = (flowDesc->bc[0][ZPLUS].type == PERIODIC)? DMDA_BOUNDARY_PERIODIC : DMDA_BOUNDARY_GHOSTED;

	// Create distributed array data structures
	// pressure
	switch(dim)
	{
		case 2: ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_STAR, mesh->nx, mesh->ny, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, &pda); CHKERRQ(ierr);
		        break;
		case 3: ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_STAR, mesh->nx, mesh->ny, mesh->nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 1, 1, NULL, NULL, NULL, &pda); CHKERRQ(ierr);
		        break;
	}
	// packed DMs
	ierr = DMCompositeCreate(PETSC_COMM_WORLD, &lambdaPack); CHKERRQ(ierr);
	ierr = DMCompositeAddDM(lambdaPack, pda); CHKERRQ(ierr);
	// velocity fluxes
	switch(dim)
	{
		case 2: ierr = DMDACreate2d(PETSC_COMM_WORLD, bx, by, DMDA_STENCIL_STAR, mesh->nx, mesh->ny, PETSC_DECIDE, PETSC_DECIDE, 2, 1, NULL, NULL, &qda); CHKERRQ(ierr);
		        break;
		case 3: ierr = DMDACreate3d(PETSC_COMM_WORLD, bx, by, bz, DMDA_STENCIL_STAR, mesh->nx, mesh->ny, mesh->nz, PETSC_DECIDE, PETSC_DECIDE, PETSC_DECIDE, 3, 1, NULL, NULL, NULL, &qda); CHKERRQ(ierr);
		        break;
	}


	return 0;
}