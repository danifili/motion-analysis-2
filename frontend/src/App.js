import React, { Component } from 'react';
import './App.css';
import BezierCurve from './BezierCurve';

class App extends Component {

  constructor(props) {
    super(props);
     this.state = { image: null };
     this.onDrop = this.onDrop.bind(this);
  }

  onImageChange = (event) => {
    if (event.target.files && event.target.files[0]) {
      this.setState({
        image: URL.createObjectURL(event.target.files[0])
      });
      console.log(this.state.image)
    }
   }

  onDrop(picture) {
      this.setState({
          pictures: this.state.pictures.concat(picture),
      });
      console.log(picture)
  }
  render() {
    return (
      <div className="App">
        <div style={{width:"40%", margin:"auto"}}>
          <BezierCurve img={this.state.image}/>
        </div>
        <input type="file" onChange={this.onImageChange}/>
      </div>
    );
  }
}

export default App;
