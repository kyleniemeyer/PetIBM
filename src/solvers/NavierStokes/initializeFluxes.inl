template <>
PetscErrorCode NavierStokesSolver<2>::initializeFluxes()
{
	PetscErrorCode ierr;
	
	if(simParams->restart)
	{
		ierr = readFluxes(); CHKERRQ(ierr);
	}
	else
	{
		PetscInt  mstart, nstart, m, n;
		PetscReal ***qArray;
		PetscReal initVel[2]  = {flowDesc->initialVelocity[0], flowDesc->initialVelocity[1]};
		PetscReal initPert[2] = {flowDesc->initialPerturbation[0], flowDesc->initialPerturbation[1]};
		PetscReal width[2]    = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0]};

		ierr = DMDAVecGetArrayDOF(qda, q, &qArray); CHKERRQ(ierr);
		ierr = DMDAGetCorners(qda, &mstart, &nstart, NULL, &m, &n, NULL); CHKERRQ(ierr);

		// Set interior values for fluxes
		for(PetscInt j=nstart; j<nstart+n; j++)
		{
			for(PetscInt i=mstart; i<mstart+m; i++)
			{
				PetscReal x, y;

				// u-fluxes
				x = -PETSC_PI + 2*PETSC_PI*(mesh->x[i+1] - mesh->x[0])/width[0];
				y = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1];
				
				qArray[j][i][0] = (initVel[0] + initPert[0]*cos(0.5*x)*sin(y))*mesh->dy[j];
			
				// v-fluxes
				x = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
				y = -PETSC_PI + 2*PETSC_PI*(mesh->y[j+1] - mesh->y[0])/width[1];
				
				qArray[j][i][1] = (initVel[1] + initPert[1]*sin(x)*cos(0.5*y))*mesh->dx[i];
			}
		}

		ierr = DMDAVecRestoreArrayDOF(qda, q, &qArray); CHKERRQ(ierr);
	}
	
	ierr = DMGlobalToLocalBegin(qda, q, INSERT_VALUES, qLocal); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(qda, q, INSERT_VALUES, qLocal); CHKERRQ(ierr);

	return 0;
}

template <>
PetscErrorCode NavierStokesSolver<3>::initializeFluxes()
{
	PetscErrorCode ierr;

	if(simParams->restart)
	{
		ierr = readFluxes(); CHKERRQ(ierr);
	}
	else
	{
		PetscInt  mstart, nstart, pstart, m, n, p;
		PetscReal ****qArray;
		PetscReal initVel[3]  = {flowDesc->initialVelocity[0], flowDesc->initialVelocity[1], flowDesc->initialVelocity[2]};
		PetscReal initPert[3] = {flowDesc->initialPerturbation[0], flowDesc->initialPerturbation[1], flowDesc->initialPerturbation[2]};
		PetscReal width[3]    = {mesh->x[mesh->nx] - mesh->x[0], mesh->y[mesh->ny] - mesh->y[0], mesh->z[mesh->nz] - mesh->z[0]};

		ierr = DMDAVecGetArrayDOF(qda, q, &qArray); CHKERRQ(ierr);
		ierr = DMDAGetCorners(qda, &mstart, &nstart, &pstart, &m, &n, &p); CHKERRQ(ierr);
		
		// Set interior values for fluxes
		for(PetscInt k=pstart; k<pstart+p; k++)
		{
			for(PetscInt j=nstart; j<nstart+n; j++)
			{
				for(PetscInt i=mstart; i<mstart+m; i++)
				{
					PetscReal x, y, z;

					// u-fluxes
					x = -PETSC_PI + 2*PETSC_PI*(mesh->x[i+1] - mesh->x[0])/width[0];
				    y = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1];
					z = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
					
					qArray[k][j][i][0] = (initVel[0] + initPert[0]*cos(0.5*x)*sin(y)*sin(z))*(mesh->dy[j]*mesh->dz[k]);
					
					// v-fluxes
					x = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
				    y = -PETSC_PI + 2*PETSC_PI*(mesh->y[j+1] - mesh->y[0])/width[1],
				    z = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->z[k]+mesh->z[k+1]) - mesh->z[0])/width[2];
				
					qArray[k][j][i][1] = (initVel[1] + initPert[1]*sin(x)*cos(0.5*y)*sin(z))*(mesh->dx[i]*mesh->dz[k]);

					// w-fluxes					
					x = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->x[i]+mesh->x[i+1]) - mesh->x[0])/width[0],
					y = -PETSC_PI + 2*PETSC_PI*(0.5*(mesh->y[j]+mesh->y[j+1]) - mesh->y[0])/width[1],
					z = -PETSC_PI + 2*PETSC_PI*(mesh->z[k+1] - mesh->z[0])/width[2];
					
					qArray[k][j][i][2] = (initVel[2] + initPert[2]*sin(x)*sin(y)*cos(0.5*z))*(mesh->dx[i]*mesh->dy[j]);
				}
			}
		}
		
		ierr = DMDAVecRestoreArrayDOF(qda, q, &qArray); CHKERRQ(ierr);
	}
	
	ierr = DMGlobalToLocalBegin(qda, q, INSERT_VALUES, qLocal); CHKERRQ(ierr);
	ierr = DMGlobalToLocalEnd(qda, q, INSERT_VALUES, qLocal); CHKERRQ(ierr);

	return 0;
}
