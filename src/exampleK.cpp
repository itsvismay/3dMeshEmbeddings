#include <vector>

//#include <Qt3DIncludes.h>
#include <GaussIncludes.h>
#include <ForceSpring.h>
#include <FEMIncludes.h>
#include <PhysicalSystemParticles.h>
//Any extra things I need such as constraints
#include <ConstraintFixedPoint.h>
#include <TimeStepperEulerImplicitLinear.h>
#include <TimeStepperEulerImplicit.h>
#include <AssemblerParallel.h>

#include <igl/get_seconds.h>
#include <igl/writeDMAT.h>
#include <igl/viewer/Viewer.h>
#include <igl/readMESH.h>
#include <igl/unproject_onto_mesh.h>

#include <stdlib.h>
#include <iostream>
#include <string>
#include <sstream>
#include <map>

using namespace Gauss;
using namespace FEM;
using namespace ParticleSystem; //For Force Spring

/* Tetrahedral finite elements */

//typedef physical entities I need

typedef PhysicalSystemFEM<double, NeohookeanTet> NeohookeanTets;

typedef World<double,
                        std::tuple<PhysicalSystemParticleSingle<double> *, NeohookeanTets *>,
                        std::tuple<ForceSpringFEMParticle<double> *>,
                        std::tuple<ConstraintFixedPoint<double> *> > MyWorld;
typedef TimeStepperEulerImplicit<double, AssemblerEigenSparseMatrix<double>,
 AssemblerEigenVector<double> > MyTimeStepper;

// Mesh
Eigen::MatrixXd V; // Verts
Eigen::MatrixXi T; // Tet indices
Eigen::MatrixXi F; // Face indices

Eigen::MatrixXd V1; // Verts
Eigen::MatrixXi T1; // Tet indices
Eigen::MatrixXi F1; // Face indices

Eigen::MatrixXd N;

// Mouse/Viewer state
Eigen::RowVector3f last_mouse;
Eigen::RowVector3d dragged_pos;
bool is_dragging = false;
int dragged_vert = 0;
int current_frame = 0;


Eigen::MatrixXd createEmbeddedMesh(Eigen::MatrixXd EmbeddingV, NeohookeanTets *tets, MyWorld &world) {
    Eigen::MatrixXd Nc;
    Nc.resize(3*EmbeddingV.rows(), 3*tets->getImpl().getV().rows()); //Embedded x originalMesh

    Nc.setZero();
    std::cout<<Nc.rows()<<", "<<Nc.cols()<<std::endl<<std::endl;
    // std::cout<<EmbeddingV<<std::endl<<std::endl;

    Eigen::VectorXd reachedThisShapeFunction(3*tets->getImpl().getV().rows());

    for(unsigned int i=0; i<EmbeddingV.rows(); ++i)
    {
        // std::cout<<"Row i "<<i<<std::endl;
        // std::cout<<EmbeddingV.row(i)<<std::endl;
        reachedThisShapeFunction.setZero();

        for(auto element: tets->getImpl().getElements())
        {
            int i0 = (element->EnergyKineticNonLumped::m_qDofs)[0]->getLocalId();
            int i1 = (element->EnergyKineticNonLumped::m_qDofs)[1]->getLocalId();
            int i2 = (element->EnergyKineticNonLumped::m_qDofs)[2]->getLocalId();
            int i3 = (element->EnergyKineticNonLumped::m_qDofs)[3]->getLocalId();

            double ro_array[3];
            Eigen::RowVectorXd::Map(ro_array, EmbeddingV.row(i).cols()) = EmbeddingV.row(i);
            double phi0 = element->EnergyKineticNonLumped::template phi<0>(ro_array);
            double phi1 = element->EnergyKineticNonLumped::template phi<1>(ro_array);
            double phi2 = element->EnergyKineticNonLumped::template phi<2>(ro_array);
            double phi3 = element->EnergyKineticNonLumped::template phi<3>(ro_array);
            if((phi0<-1e-6 || phi0>(1+1e-6)) || (phi1<-1e-6 || phi1>(1+1e-6)) || (phi2<-1e-6 || phi2>(1+1e-6)) || (phi3<-1e-6 || phi3>(1+1e-6)))
            {
                // std::cout<<"phi val: "<<phi0<<", "<<phi1<<", "<<phi2<<", "<<phi3<<std::endl;
                phi0 = 0;
                phi1 = 0;
                phi2 = 0;
                phi3 = 0;
            }


            // std::cout<<"phi val: "<<phi0<<", "<<phi1<<", "<<phi2<<", "<<phi3<<std::endl;
            // std::cout<<"local Ids: "<<i0<<", "<<i1<<", "<<i2<<", "<<i3<<std::endl;
            // std::cout<<"DOFs: "<<
            // mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[0], world.getState()).transpose()<<", "<<
            // mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[1], world.getState()).transpose()<<", "<<
            // mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[2], world.getState()).transpose()<<", "<<
            // mapDOFEigen(*(element->EnergyKineticNonLumped::m_qDofs)[3], world.getState()).transpose()<<
            // std::endl;

            // std::cout<<reachedThisShapeFunction.transpose()<<std::endl;
            if(reachedThisShapeFunction(i0)<1e-6){
                // std::cout<<"i0 "<<i0<<std::endl;
                Nc(3*i,   i0 )   += phi0;
                Nc(3*i+1, i0+1 ) += phi0;
                Nc(3*i+2, i0+2 ) += phi0;
                reachedThisShapeFunction(i0) = phi0;
            }

            if(reachedThisShapeFunction(i1)<1e-6){
                // std::cout<<"i1 "<<i1<<std::endl;
                Nc(3*i,   i1 )   += phi1;
                Nc(3*i+1, i1+1 ) += phi1;
                Nc(3*i+2, i1+2 ) += phi1; 
                reachedThisShapeFunction(i1) = phi1;
            }

            if(reachedThisShapeFunction(i2)<1e-6){
                // std::cout<<"i2 "<<i2<<std::endl;
                Nc(3*i,   i2 )   += phi2;
                Nc(3*i+1, i2+1 ) += phi2;
                Nc(3*i+2, i2+2 ) += phi2; 
                reachedThisShapeFunction(i2) = phi2;
            }

            if(reachedThisShapeFunction(i3)<1e-6){
                // std::cout<<"i3 "<<i3<<std::endl;
                Nc(3*i,   i3 )   += phi3;
                Nc(3*i+1, i3+1 ) += phi3;
                Nc(3*i+2, i3+2 ) += phi3;  
                reachedThisShapeFunction(i3) = phi3;
            } 
        } 
        // std::cout<<"---------------------"<<std::endl;
        // std::cout<<std::endl;    
    }

    // std::cout<<Nc<<std::endl;
    return Nc;
}

// Todo put this in utilities
Eigen::MatrixXd getCurrentVertPositions(MyWorld &world, NeohookeanTets *tets) {
    // Eigen::Map<Eigen::MatrixXd> q(mapStateEigen<0>(world).data(), V.cols(), V.rows()); // Get displacements only
    auto q = mapDOFEigen(tets->getQ(), world);
    Eigen::VectorXd q1 = N*q;
    Eigen::Map<Eigen::MatrixXd> dV(q1.data(), V1.cols(), V1.rows()); // Get displacements only

    return V1 + dV.transpose();
}


void reset_world (MyWorld &world) {
        auto q = mapStateEigen(world); // TODO is this necessary?
        q.setZero();
}

int main(int argc, char **argv) {

    //Setup Physics
    MyWorld world;

    igl::readMESH(argv[1], V, T, F);
    std::cout<<"Test Neohookean FEM \n";
    NeohookeanTets *tets = new NeohookeanTets(V,T);
    Eigen::Vector3d x(0,0,0);
    for(auto element: tets->getImpl().getElements()) {
        element->EnergyKineticNonLumped::setDensity(100.0);//1000.0);
        element->setParameters(1000000, 0.45);
    }

    // // Pinned particle to attach spring for dragging
    PhysicalSystemParticleSingle<double> *pinned_point = new PhysicalSystemParticleSingle<double>();
    pinned_point->getImpl().setMass(100000000); //10000000
    auto fem_attached_pos = PosFEM<double>(&tets->getQ()[0],0, &tets->getImpl().getV());

    double spring_stiffness = 150.0;
    double spring_rest_length = 0.1;
    ForceSpringFEMParticle<double> *forceSpring = new ForceSpringFEMParticle<double>(fem_attached_pos, // TODO compare getV to V. Get rid of double use of index
                                                                                     PosParticle<double>(&pinned_point->getQ()),
                                                                                     spring_rest_length, 0.0);

    world.addSystem(pinned_point);
    world.addForce(forceSpring);
    world.addSystem(tets);
    fixDisplacementMin(world, tets, 0);
    world.finalize(); //After this all we're ready to go (clean up the interface a bit later)

    reset_world(world);

    MyTimeStepper stepper(0.05, 10);
    stepper.step(world);

    igl::readMESH(argv[2], V1, T1, F1);
    N = createEmbeddedMesh(V1, tets, world);


    /** libigl display stuff **/
    igl::viewer::Viewer viewer;
    double tot_time = 0.0;
    viewer.callback_pre_draw = [&](igl::viewer::Viewer & viewer)
    {
        // predefined colors
        const Eigen::RowVector3d orange(1.0,0.7,0.2);
        const Eigen::RowVector3d yellow(1.0,0.9,0.2);
        const Eigen::RowVector3d blue(0.2,0.3,0.8);
        const Eigen::RowVector3d green(0.2,0.6,0.3);

        if(is_dragging) {
            auto q_part = mapDOFEigen(pinned_point->getQ(), world);
            Eigen::MatrixXd part_pos = q_part;

            viewer.data.set_points(part_pos.transpose(), orange);
        } else {
            Eigen::MatrixXd part_pos;
            viewer.data.set_points(part_pos, orange);
        }

        if(viewer.core.is_animating)
        {
            // Save Current configuration
            double start = igl::get_seconds();
            stepper.step(world);
            tot_time += igl::get_seconds() - start;

            Eigen::MatrixXd newV = getCurrentVertPositions(world, tets);
            viewer.data.set_vertices(newV);
            viewer.data.compute_normals();

            current_frame++;
            // std::cout << "Average timestep: " << (tot_time / current_frame) << std::endl;
        }
        return false;
    };

    viewer.callback_key_pressed = [&](igl::viewer::Viewer &, unsigned int key, int mod)
    {
        switch(key)
        {
            case 'P':
            case 'p':
            {
                viewer.core.is_animating = !viewer.core.is_animating;
                break;
            }
            case 'r':
            case 'R':
            {
                reset_world(world);
                break;
            }
            default:
            return false;
        }
        return true;
    };

    viewer.callback_mouse_down = [&](igl::viewer::Viewer&, int, int)->bool
    {
        Eigen::MatrixXd curV = getCurrentVertPositions(world, tets);
        last_mouse = Eigen::RowVector3f(viewer.current_mouse_x,viewer.core.viewport(3)-viewer.current_mouse_y,0);

        // Find closest point on mesh to mouse position
        int fid;
        Eigen::Vector3f bary;
        if(igl::unproject_onto_mesh(
            last_mouse.head(2),
            viewer.core.view * viewer.core.model,
            viewer.core.proj,
            viewer.core.viewport,
            curV, F,
            fid, bary))
        {
            long c;
            bary.maxCoeff(&c);
            dragged_pos = curV.row(F(fid,c)) + Eigen::RowVector3d(0.001,0.0,0.0); //Epsilon offset so we don't div by 0
            dragged_vert = F(fid,c);

            // Update the system
            is_dragging = true;
            forceSpring->getImpl().setStiffness(spring_stiffness);
            auto pinned_q = mapDOFEigen(pinned_point->getQ(), world);
            pinned_q = dragged_pos;//(dragged_pos).cast<double>(); // necessary?
            std::cout << "Dragging vert: " << dragged_pos << std::endl;

            fem_attached_pos = PosFEM<double>(&tets->getQ()[dragged_vert],dragged_vert, &tets->getImpl().getV());
            forceSpring->getImpl().setPosition0(fem_attached_pos);

            return true;
        }

        return false; // TODO false vs true??
    };

    viewer.callback_mouse_up = [&](igl::viewer::Viewer&, int, int)->bool
    {
        is_dragging = false;
        forceSpring->getImpl().setStiffness(0.0);
        // auto pinned_q = mapDOFEigen(pinned_point->getQ(), world);
        // pinned_q = Eigen::RowVector3d(-10000.0,0.0,0.0); // Move the point out of view

        return false;
    };

    viewer.callback_mouse_move = [&](igl::viewer::Viewer &, int,int)->bool
    {
        if(is_dragging) {
            Eigen::RowVector3f drag_mouse(
                viewer.current_mouse_x,
                viewer.core.viewport(3) - viewer.current_mouse_y,
                last_mouse(2));

            Eigen::RowVector3f drag_scene,last_scene;

            igl::unproject(
                drag_mouse,
                viewer.core.view*viewer.core.model,
                viewer.core.proj,
                viewer.core.viewport,
                drag_scene);
            igl::unproject(
                last_mouse,
                viewer.core.view*viewer.core.model,
                viewer.core.proj,
                viewer.core.viewport,
                last_scene);

            dragged_pos += ((drag_scene-last_scene)*4.5).cast<double>(); //TODO why do I need to fine tune this
            // dragged_pos += (drag_mouse-last_mouse).cast<double>();
            last_mouse = drag_mouse;

            // Update the system
            // TODO dedupe this
            auto pinned_q = mapDOFEigen(pinned_point->getQ(), world);
            pinned_q = dragged_pos; //(dragged_pos).cast<double>(); // necessary?
        }

        return false;
    };

    viewer.data.set_mesh(V1,F1);
    viewer.core.show_lines = true;
    viewer.core.invert_normals = true;
    viewer.core.is_animating = false;
    viewer.data.face_based = true;

    viewer.launch();
    return EXIT_SUCCESS;
}