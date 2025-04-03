surfaceScalarField phiv("phiv", 0.0 * phiv_pos);

surfaceVectorField phiU("phiU", 0.0 * aphiv_pos * rhoU_pos);
phiU.setOriented(true);
            
surfaceVectorField phiUp("phiUp", 0.0 * p_pos * mesh.Sf());
phiUp.setOriented(true);
        
surfaceScalarField phiEp("phiEp", 0.0 * phiv_pos * p_pos);
phiEp.setOriented(false);

