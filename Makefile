# project name (generate executable with this name)
TARGET_PENDULE   = optimal_pendule
TARGET_NONHOLO   = optimal_nonholo
TARGET_NONHOLO2   = optimal_nonholo2
TARGET_LEG2DTHETA   = optimal_leg2DTheta
TARGET_LEG2DH   = optimal_leg2DH
# TARGET_LEG2D   = optimal_leg2D
TARGET_LIPM   = optimal_lipm
TARGET_LIPM2   = optimal_lipm2


CC       = g++
#CC       = ccache g++

# compiling flags here
FLAGS = -g -O3 #-D_PROBACTIONTHRES #-D_NOSERIAL #-D_MEMOPTIM

# change these to set the proper directories where each files shoould be
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin
INCLDIR  = include

#LIBS= -lm -lboost_graph -lboost_serialization
LIBS= -g -lm -L/usr/local/lib -lboost_graph -lboost_serialization

#SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
SOURCES = $(SRCDIR)/Controller.cpp $(SRCDIR)/Learner.cpp $(SRCDIR)/Manager.cpp $(SRCDIR)/Optimizer.cpp $(SRCDIR)/Specification.cpp
INCLUDES := $(wildcard $(INCLDIR)/*.h)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
rm       = rm -f

all: $(BINDIR)/$(TARGET_PENDULE) $(BINDIR)/$(TARGET_NONHOLO) $(BINDIR)/$(TARGET_NONHOLO2) $(BINDIR)/$(TARGET_LEG2DTHETA) $(BINDIR)/$(TARGET_LEG2DH) $(BINDIR)/$(TARGET_LIPM) $(BINDIR)/$(TARGET_LIPM2)

$(BINDIR)/$(TARGET_PENDULE): $(OBJECTS) $(OBJDIR)/PendulumDescription.o $(OBJDIR)/$(TARGET_PENDULE).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/PendulumDescription.o $(OBJDIR)/$(TARGET_PENDULE).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_PENDULE).o: $(SRCDIR)/$(TARGET_PENDULE).cpp  $(INCLDIR)/Models/Pendulum/PendulumDescription.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/PendulumDescription.o: $(INCLDIR)/Models/Pendulum/PendulumDescription.cpp  $(INCLDIR)/Models/Pendulum/PendulumDescription.h  $(INCLDIR)/Models/Pendulum/TodorovConfig.h  $(INCLDIR)/Models/Pendulum/WalkingConfig.h  $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)


$(BINDIR)/$(TARGET_NONHOLO): $(OBJECTS) $(OBJDIR)/NonHoloDescription.o $(OBJDIR)/$(TARGET_NONHOLO).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/NonHoloDescription.o $(OBJDIR)/$(TARGET_NONHOLO).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_NONHOLO).o: $(SRCDIR)/$(TARGET_NONHOLO).cpp  $(INCLDIR)/Models/NonHolo/NonHoloDescription.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/NonHoloDescription.o: $(INCLDIR)/Models/NonHolo/NonHoloDescription.cpp  $(INCLDIR)/Models/NonHolo/NonHoloDescription.h  $(INCLDIR)/Models/NonHolo/PaperConfig.h $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)


$(BINDIR)/$(TARGET_NONHOLO2): $(OBJECTS) $(OBJDIR)/NonHoloDescription2.o $(OBJDIR)/$(TARGET_NONHOLO2).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/NonHoloDescription2.o $(OBJDIR)/$(TARGET_NONHOLO2).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_NONHOLO2).o: $(SRCDIR)/$(TARGET_NONHOLO2).cpp  $(INCLDIR)/Models/NonHolo2/NonHoloDescription2.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/NonHoloDescription2.o: $(INCLDIR)/Models/NonHolo2/NonHoloDescription2.cpp  $(INCLDIR)/Models/NonHolo2/NonHoloDescription2.h  $(INCLDIR)/Models/NonHolo2/PaperConfig2.h $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)




$(BINDIR)/$(TARGET_LEG2DTHETA): $(OBJECTS) $(OBJDIR)/Leg2DThetaDescription.o $(OBJDIR)/$(TARGET_LEG2DTHETA).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/Leg2DThetaDescription.o $(OBJDIR)/$(TARGET_LEG2DTHETA).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_LEG2DTHETA).o: $(SRCDIR)/$(TARGET_LEG2DTHETA).cpp  $(INCLDIR)/Models/Leg2DTheta/Leg2DThetaDescription.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Leg2DThetaDescription.o: $(INCLDIR)/Models/Leg2DTheta/Leg2DThetaDescription.cpp  $(INCLDIR)/Models/Leg2DTheta/Leg2DThetaDescription.h  $(INCLDIR)/Models/Leg2DTheta/Config.h $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)


$(BINDIR)/$(TARGET_LEG2DH): $(OBJECTS) $(OBJDIR)/Leg2DHDescription.o $(OBJDIR)/$(TARGET_LEG2DH).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/Leg2DHDescription.o $(OBJDIR)/$(TARGET_LEG2DH).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_LEG2DH).o: $(SRCDIR)/$(TARGET_LEG2DH).cpp  $(INCLDIR)/Models/Leg2DH/Leg2DHDescription.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Leg2DHDescription.o: $(INCLDIR)/Models/Leg2DH/Leg2DHDescription.cpp  $(INCLDIR)/Models/Leg2DH/Leg2DHDescription.h  $(INCLDIR)/Models/Leg2DH/Config.h $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)




$(BINDIR)/$(TARGET_LIPM): $(OBJECTS) $(OBJDIR)/LipmDescription.o $(OBJDIR)/$(TARGET_LIPM).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/LipmDescription.o $(OBJDIR)/$(TARGET_LIPM).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_LIPM).o: $(SRCDIR)/$(TARGET_LIPM).cpp  $(INCLDIR)/Models/Lipm/LipmDescription.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/LipmDescription.o: $(INCLDIR)/Models/Lipm/LipmDescription.cpp  $(INCLDIR)/Models/Lipm/LipmDescription.h  $(INCLDIR)/Models/Lipm/LipmConfig.h $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(BINDIR)/$(TARGET_LIPM2): $(OBJECTS) $(OBJDIR)/Lipm2Description.o $(OBJDIR)/$(TARGET_LIPM2).o
	$(CC) -o $@ $(OBJECTS) $(OBJDIR)/Lipm2Description.o $(OBJDIR)/$(TARGET_LIPM2).o $(FLAGS) $(LIBS)

$(OBJDIR)/$(TARGET_LIPM2).o: $(SRCDIR)/$(TARGET_LIPM2).cpp  $(INCLDIR)/Models/Lipm2/Lipm2Description.h 
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Lipm2Description.o: $(INCLDIR)/Models/Lipm2/Lipm2Description.cpp  $(INCLDIR)/Models/Lipm2/Lipm2Description.h  $(INCLDIR)/Models/Lipm2/Lipm2Config.h $(INCLDIR)/Manager.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)



$(OBJDIR)/Manager.o: $(SRCDIR)/Manager.cpp  $(INCLDIR)/Manager.h  $(INCLDIR)/Controller.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Controller.o:  $(SRCDIR)/Controller.cpp  $(INCLDIR)/Controller.h  $(INCLDIR)/Learner.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Learner.o: $(SRCDIR)/Learner.cpp  $(INCLDIR)/Learner.h $(INCLDIR)/Optimizer.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Optimizer.o:  $(SRCDIR)/Optimizer.cpp  $(INCLDIR)/Optimizer.h  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)

$(OBJDIR)/Specification.o: $(SRCDIR)/Specification.cpp  $(INCLDIR)/Specification.h
	$(CC) -I$(INCLDIR) $< -o $@ -c $(FLAGS) $(LIBS)



#$(OBJECTS): $(SOURCES)
#	$(CXX) -I$(INCLDIR) $? -c $(FLAGS) $(LIBS)

#$(SOURCES): $(INCLUDES)

.PHONEY: clean

clean:
	@$(rm) $(OBJECTS)

.PHONEY: remove
remove: clean
	@$(rm) $(BINDIR)/$(TARGET_PENDULE)
