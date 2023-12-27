/**
 * CS 114 final project
 * fluid simulation based on smoothed particles hydrodynamics with grid structures
 * Yuxuan Zhu
 */

class SPHFluid {
  //keep contructor code
  constructor() {
    // Physical attrs
    this.numParticles = 800;
    this.viscousity = 900 * 5;
    this.particleMass = 500 * .13;
    this.h = 16;
    this.stiffness = 400 * 5;
    this.gravityConst = 120000 * 9.82;
    this.dt = 0.0004;

    //grid attrs
    this.gridSize = this.h; //?=h
    this.grid = {};

    this.particles_ = [];
    this.particlePositions_ = [this.numParticles];
    this.fireParticles_ = false;
    
    this.initParticles();
  }

  get particles() { return this.particles_; }
  get particlePositions() { return this.particlePositions_; }

  fireParticles() {
    this.fireParticles_ = !this.fireParticles_;
  }

  //helper function for adding grid
  getgridkey(position) {
    const x = Math.floor(position.x / this.gridSize);
    const y = Math.floor(position.y / this.gridSize);
    return `${x},${y}`;
  }

  //keep initParticles code
  initParticles() {
    this.particles_ = [];

    const geometry = new THREE.CircleGeometry(WaterEngine.PARTICLE_RADIUS, 16);
    const material = new THREE.MeshBasicMaterial({color: 0x2FA1D6});
    // Set starting positions
    let k = 0;
    let j = 0;

    for (let i = 0; i < this.numParticles; i++) {
      this.particles_.push({
        position: new THREE.Vector3(0, 0, 0),
        vel: new THREE.Vector3(0, 0, 0),
        pressure: 0,
        density: 0,
        viscousityForce: new THREE.Vector3(0, 0, 0),
        pressureForce: new THREE.Vector3(0, 0, 0),
        gravityForce: new THREE.Vector3(0, 0, 0),
        otherForce: new THREE.Vector3(0, 0, 0),
      });

      if (i % 40 === 0) {
        k++;
        j = 0;
      }
      j++;

      this.particles_[i].position.set(
          WaterEngine.START_OFFSET_X + j * this.h / 2, WaterEngine.START_OFFSET_Y + k * this.h / 2, 0);
      this.particlePositions_[i] = this.particles_[i].position;

      //add particles to the grid
      const gridkey = this.getgridkey(this.particles_[i].position);
      if (!this.grid[gridkey]) {
        this.grid[gridkey] = [];
      }
      this.grid[gridkey].push(i);
    }
  }

  //rewrite
  calculateDensityAndPressure() {
    const h2 = Math.pow(this.h, 2.0);
    const h9 = Math.pow(this.h, 9.0);
    const pre_poly6 = this.particleMass * 315 / (64 * Math.PI * h9)
    const density = 998 //parameter from the original code
    for (let i = 0; i < this.particles_.length; i++) {
      let densitySum = 0;

      //get the neighboring cells
      const gridkey = this.getgridkey(this.particles_[i].position);
      const [x, y] = gridkey.split(",").map(Number);
      const neighbors = [];
      for (let dx = -1; dx <= 1; dx++) {
          for (let dy = -1; dy <= 1; dy++) {
              const neighborKey = `${x+dx},${y+dy}`;
              if (this.grid[neighborKey]) {
                  neighbors.push(...this.grid[neighborKey]);
              }
          }
      }
      //only find neighbor 
      for (let j = 0 ; j < neighbors.length; j++) {
        let diffVec = new THREE.Vector3(0, 0, 0);
        diffVec.subVectors(this.particles_[i].position, this.particles_[neighbors[j]].position);
        let lengVec = diffVec.length();

        if (lengVec <= this.h) {
          let lengVec2 = Math.pow(lengVec, 2); // lengVec^2
          let diff = Math.pow((h2 - lengVec2), 3.0);

          densitySum += pre_poly6 * diff; // Sum up density
        }
      }

      this.particles_[i].density = densitySum;
      this.particles_[i].pressure = this.stiffness * (densitySum - density);

      //flip gravity
      if (this.fireParticles_) {
        this.particles_[i].otherForce = new THREE.Vector3(0, 50000, 0);
      } else {
        this.particles_[i].otherForce = new THREE.Vector3(0, 0, 0);
      }
    }
  }

  //rewrite
  calculateForces() {
    const h6 = Math.pow(this.h, 6.0);
    const pre_spiky = - this.particleMass * 45 / (Math.PI * h6)
    const pre_lapla = this.viscousity * this.particleMass * 45 / (Math.PI * h6)

    for (let i = 0; i < this.numParticles; i++) {
      let gravity = new THREE.Vector3(0, -this.gravityConst * this.particles_[i].density, 0);
      let pressure = new THREE.Vector3(0, 0, 0);
      let viscousity = new THREE.Vector3(0, 0, 0);

      const gridkey = this.getgridkey(this.particles_[i].position);
      const [x, y] = gridkey.split(",").map(Number);
      const neighbors = [];
      for (let dx = -1; dx <= 1; dx++) {
          for (let dy = -1; dy <= 1; dy++) {
              const neighborKey = `${x+dx},${y+dy}`;
              if (this.grid[neighborKey]) {
                  neighbors.push(...this.grid[neighborKey]);
              }
          }
      }
      //only find neighbor 
      for (let j = 0; j < neighbors.length; j++) {
        if (i === j) {
          continue;
        }
        let diffVec = new THREE.Vector3(0, 0, 0);
        diffVec.subVectors(this.particles_[i].position, this.particles_[neighbors[j]].position);
        let lengVec = diffVec.length();
        if (lengVec <= this.h) {
          let diff = this.h - lengVec;
          let pre_average = (this.particles_[i].pressure + this.particles_[neighbors[j]].pressure) / 2.0;
          let pre_temp = pre_spiky * Math.pow(diff, 3.0) * pre_average / lengVec / this.particles_[neighbors[j]].density;
          let p = new THREE.Vector3(pre_temp * diffVec.x, pre_temp * diffVec.y, 0);
          pressure.add(p);

          let diffVel = new THREE.Vector3(0, 0, 0);
          diffVel.subVectors(this.particles_[neighbors[j]].vel, this.particles_[i].vel);
          let v = diffVel.multiplyScalar(pre_lapla * diff / this.particles_[neighbors[j]].density)
          viscousity.add(v);
        }
      }
      this.particles_[i].gravityForce.set(gravity.x, gravity.y, gravity.z);
      this.particles_[i].viscousityForce.set(viscousity.x, viscousity.y, viscousity.z);
      this.particles_[i].pressureForce.set(pressure.x, pressure.y, pressure.z);
    }
  }

  calculateAcceleration() {
    this.calculateDensityAndPressure();
    this.calculateForces();
  }

  //rewrite
  idle() {
    let newPos = new THREE.Vector3(0, 0, 0);
    let newVel = new THREE.Vector3(0, 0, 0);
    let newPositions = [];

    for (let i = 0; i < this.particles_.length; i++) {
      newPos.addVectors(this.particles_[i].gravityForce, this.particles_[i].viscousityForce);
      newPos.add(this.particles_[i].pressureForce);
      newPos.add(this.particles_[i].otherForce)
          newPos.multiplyScalar((this.dt * this.dt) / (2 * this.particles_[i].density));
      newPos.add(this.particles_[i].vel.multiplyScalar(this.dt));
      newPos.add(this.particles_[i].position);

      newVel.subVectors(newPos, this.particles_[i].position);
      newVel.multiplyScalar(1 / this.dt);

      this.particles_[i].position.set(newPos.x, newPos.y, newPos.z);
      this.particles_[i].vel.set(newVel.x, newVel.y, newVel.z);
      newPositions.push(this.particles_[i].position);
      this.checkBoundaries(this.particles_[i]);
    }
    this.particlePositions_ = newPositions;

    //clear and refill grid
    this.grid = {};
    for (let i = 0; i < this.particles_.length; i++) {
      const gridkey = this.getgridkey(this.particles_[i].position);
      if (!this.grid[gridkey]) {
        this.grid[gridkey] = [];
      }
    this.grid[gridkey].push(i);
    }
    // console.log(this.grid);
  }

  //update bounding (parameter update)
  checkBoundaries(particle) {
    if (particle.position.x < WaterEngine.PARTICLE_RADIUS + 1) {
      particle.vel.x = -0.8 * particle.vel.x;
      particle.position.x = WaterEngine.PARTICLE_RADIUS + 2;
    }

    else if (particle.position.x > WaterEngine.SQUARE_SIZE - WaterEngine.PARTICLE_RADIUS - 1) {
      particle.vel.x = -0.8 * particle.vel.x;
      particle.position.x = WaterEngine.SQUARE_SIZE - WaterEngine.PARTICLE_RADIUS - 2;
    }

    if (particle.position.y < WaterEngine.PARTICLE_RADIUS + 1) {
      particle.vel.y = -0.8 * particle.vel.y;
      particle.position.y = WaterEngine.PARTICLE_RADIUS + 2;

    } else if (
      particle.position.y > WaterEngine.SQUARE_SIZE - WaterEngine.PARTICLE_RADIUS - 1) {
      particle.vel.y = -0.8 * particle.vel.y;
      particle.position.y = WaterEngine.SQUARE_SIZE - WaterEngine.PARTICLE_RADIUS - 2;
    }
  }
}
