import Kidneys from './svg/kidneys.svg';

interface Icon {
  size: number;
}
function KidneyIcon({ size }: Icon) {
  return (
    <div style={{ width: size }}>
      <img src={Kidneys} alt="kidney icon" />
    </div>
  );
}

export default KidneyIcon;
